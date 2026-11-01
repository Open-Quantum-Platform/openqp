! Fixed-reference LDA/GGA response data. Missing entries are recomputed, never
! treated as zero. The owner releases this cache at the end of its solve.
module mod_dft_gridint_response_cache
  use precision, only: fp, i8b
  implicit none
  private
  public :: response_cache_t, response_slice_t, cache_hash

  type response_slice_t
    logical :: ao_ready=.false., ao_skip=.false., dense=.true.
    logical :: kernel_ready=.false., kernel_skip=.false.
    integer :: npts=0, naos=0
    integer, allocatable :: indices(:)
    real(fp), allocatable :: ao(:), weights(:), kernel(:,:)
    real(fp) :: stats(5)=0.0_fp
  end type

  type response_cache_t
    type(response_slice_t), allocatable :: slices(:)
    integer(i8b) :: key=0, ao_key=0, max_bytes=-1, configured_limit=0, bytes=0, ao_bytes=0, kernel_bytes=0
    integer :: kernel_values=32
    integer(i8b) :: ao_limit=0, kernel_limit=0
    logical :: enabled=.false.
    integer :: ao_hits=0, kernel_hits=0
  contains
    procedure :: budget
    procedure :: prepare
    procedure :: store_ao
    procedure :: store_kernel
    procedure :: free
  end type

  interface cache_hash
    module procedure hash_reals, hash_matrix, hash_ints
  end interface
contains
  ! Hash lengths as well as values so different array partitions cannot alias.
  subroutine hash_ints(h, a)
    integer(i8b), intent(inout) :: h
    integer, intent(in) :: a(:)
    integer :: i
    h=ieor(h,int(size(a),i8b))*1099511628211_i8b
    do i=1,size(a)
      h=ieor(h,int(a(i),i8b))*1099511628211_i8b
    end do
  end subroutine
  subroutine hash_reals(h, a)
    integer(i8b), intent(inout) :: h
    real(fp), intent(in) :: a(:)
    integer :: i
    h=ieor(h,int(size(a),i8b))*1099511628211_i8b
    do i=1,size(a)
      h=ieor(h,transfer(a(i),0_i8b))*1099511628211_i8b
    end do
  end subroutine
  subroutine hash_matrix(h, a)
    integer(i8b), intent(inout) :: h
    real(fp), intent(in) :: a(:,:)
    integer :: j
    call hash_ints(h,shape(a))
    do j=1,size(a,2)
      call hash_reals(h,a(:,j))
    end do
  end subroutine

  function budget(self) result(limit)
    class(response_cache_t), intent(in) :: self
    integer(i8b) :: limit
    integer :: stat,ln
    real(fp) :: mb
    character(64) :: value
    limit=self%max_bytes
    if(limit>=0) return
    mb=256.0_fp
    call get_environment_variable('OQP_XC_RESPONSE_CACHE_MB',value,length=ln,status=stat)
    if(stat==0.and.ln>0) then
      read(value,*,iostat=stat) mb
      if(stat/=0) mb=0.0_fp
    end if
    if(.not.(mb>=0.0_fp.and.mb<=1048576.0_fp)) mb=0.0_fp
    limit=int(mb*1048576.0_fp,i8b)
  end function

  subroutine prepare(self, key, nslices, npoints, ao_key, kernel_values)
    class(response_cache_t), intent(inout) :: self
    integer(i8b), intent(in) :: key
    integer, intent(in) :: nslices,npoints
    integer(i8b), intent(in), optional :: ao_key
    integer, intent(in), optional :: kernel_values
    integer(i8b) :: limit, metadata, available, geometry_key
    integer :: stat, nv, i
    type(response_slice_t) :: sample
    limit=self%budget()
    geometry_key=key
    if(present(ao_key)) geometry_key=ao_key
    nv=32
    if(present(kernel_values)) nv=kernel_values
    self%ao_hits=0; self%kernel_hits=0
    if(allocated(self%slices)) then
      if(self%ao_key==geometry_key.and.size(self%slices)==nslices.and. &
         self%configured_limit==limit.and.self%kernel_values==nv) then
        if(self%key/=key) then
          do i=1,nslices
            if(allocated(self%slices(i)%kernel)) deallocate(self%slices(i)%kernel)
            self%slices(i)%kernel_ready=.false.
            self%slices(i)%kernel_skip=.false.
            self%slices(i)%stats=0.0_fp
          end do
          self%bytes=self%bytes-self%kernel_bytes
          self%kernel_bytes=0
          self%key=key
        end if
        return
      end if
    end if
    call self%free()
    metadata=int(storage_size(sample)/8,i8b)*int(nslices,i8b)
    if(nslices<=0.or.limit<=metadata) return
    allocate(self%slices(nslices),stat=stat)
    if(stat/=0) return
    self%key=key; self%bytes=metadata; self%configured_limit=limit
    self%ao_key=geometry_key; self%kernel_values=nv
    available=limit-metadata
    ! Reserve compact reference data before admitting the much larger AO data.
    self%kernel_limit=min(available,8_i8b*nv*int(npoints,i8b))
    self%ao_limit=available-self%kernel_limit
    self%enabled=.true.
  end subroutine

  subroutine store_ao(self,islice,skip,np,na,dense,nvec,indices,ao,weights)
    class(response_cache_t), intent(inout) :: self
    integer,intent(in) :: islice,np,na,nvec,indices(:)
    logical,intent(in) :: skip,dense
    real(fp),intent(in) :: ao(:),weights(:)
    integer(i8b) :: need,n
    integer :: stat
    if(.not.self%enabled) return
    associate(s=>self%slices(islice))
      if(s%ao_ready) return
      s%ao_skip=skip; s%npts=np; s%naos=na; s%dense=dense
      if(skip) then
        s%ao_ready=.true.
        return
      end if
      n=int(na,i8b)*np*nvec
      need=8_i8b*(n+np)
      if(.not.dense) need=need+int(storage_size(na)/8,i8b)*na
!$omp critical(oqp_response_cache_allocation)
      if(need<=self%ao_limit-self%ao_bytes) then
        allocate(s%ao(n),s%weights(np),stat=stat)
        if(stat==0.and..not.dense) allocate(s%indices(na),stat=stat)
        if(stat==0) then
          s%ao=ao(:n); s%weights=weights(:np)
          if(.not.dense) s%indices=indices(:na)
          self%ao_bytes=self%ao_bytes+need; self%bytes=self%bytes+need
          s%ao_ready=.true.
        else
          if(allocated(s%ao)) deallocate(s%ao)
          if(allocated(s%weights)) deallocate(s%weights)
          if(allocated(s%indices)) deallocate(s%indices)
        end if
      end if
!$omp end critical(oqp_response_cache_allocation)
    end associate
  end subroutine

  subroutine store_kernel(self,islice,skip,data,stats)
    class(response_cache_t), intent(inout) :: self
    integer,intent(in) :: islice
    logical,intent(in) :: skip
    real(fp),intent(in) :: data(:,:),stats(5)
    integer(i8b) :: need
    integer :: stat
    if(.not.self%enabled) return
    associate(s=>self%slices(islice))
      if(s%kernel_ready) return
      s%kernel_skip=skip; s%stats=stats
      if(skip) then
        s%kernel_ready=.true.
        return
      end if
      need=8_i8b*size(data,kind=i8b)
!$omp critical(oqp_response_cache_allocation)
      if(need<=self%kernel_limit-self%kernel_bytes) then
        allocate(s%kernel(size(data,1),size(data,2)),stat=stat)
        if(stat==0) then
          s%kernel=data
          self%kernel_bytes=self%kernel_bytes+need; self%bytes=self%bytes+need
          s%kernel_ready=.true.
        end if
      end if
!$omp end critical(oqp_response_cache_allocation)
    end associate
  end subroutine

  subroutine free(self)
    class(response_cache_t), intent(inout) :: self
    if(allocated(self%slices)) deallocate(self%slices)
    self%enabled=.false.; self%key=0; self%ao_key=0; self%bytes=0; self%configured_limit=0
    self%ao_bytes=0; self%kernel_bytes=0; self%ao_limit=0; self%kernel_limit=0
    self%ao_hits=0; self%kernel_hits=0
  end subroutine
end module
