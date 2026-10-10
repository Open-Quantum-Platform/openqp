! Immutable CPU contraction excerpts used only by the source-equivalence test.
subroutine int2_mrsf_data_t_update(this, buf)

    implicit none

    class(int2_mrsf_data_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf
    integer :: i, j, k, l, n
    real(kind=dp) :: val, xval, cval
    integer :: mythread, gpu_ierr

    mythread = buf%thread_id

    if (.not.this%tamm_dancoff) return

    if (gpu_backend_metc_enabled()) then
      ! Resident hot path: borrow d3/f3/ids/ints workspace pointers; copy this
      ! flush's ids/ints into the thread's resident scratch and launch.  No
      ! per-flush cudaMalloc/cudaFree; f3 stays resident and is downloaded once
      ! at parallel_stop.  thread is 0-based for the C ABI; buf%buf_size is the
      ! resident scratch capacity.
      call gpu_backend_metc_try_contract(this%nfocks, ubound(this%d3,2), &
                                         ubound(this%d3,3), this%nthreads, &
                                         this%cur_pass, buf%buf_size, this%d3, &
                                         buf%thread_id - 1, buf%ids(:,1:buf%ncur), &
                                         buf%ints(1:buf%ncur), buf%ncur, &
                                         this%scale_exchange, this%scale_coulomb, &
                                         .false., gpu_ierr)
      if (gpu_ierr == 0) then
        buf%ncur = 0
        return
      end if
      ! gpu_ierr /= 0 here means the GPU declined this accumulation BEFORE it
      ! became active (soft fallback); the host loop below then safely owns the
      ! whole accumulation.  An ACTIVE-session failure never reaches this point:
      ! gpu_backend_metc_try_contract aborts internally rather than allow a
      ! partial device f3 to be mixed with a host f3 (all-or-nothing contract).
    end if

    associate ( f3 => this%f3(:,:,:,:,mythread), &
                d3 => this%d3, &
                nf => this%nfocks &
      )

      do n = 1, buf%ncur
        i = buf%ids(1,n)
        j = buf%ids(2,n)
        k = buf%ids(3,n)
        l = buf%ids(4,n)
        val = buf%ints(n)

        xval = val * this%scale_exchange
        cval = val * this%scale_coulomb

        ! f3(nF,1:7,:,:) !> 1=ado2v, 2=ado1v, 3=adco1, 4=adco2, 5=ao21v, 6=aco12, 7=agdlr
        ! d3(nF,1:7,:,:) !> 1= bo2v, 2= bo1v, 3= bco1, 4= bco2, 5= o21v, 6= co12, 7= ball
        if (this%cur_pass==1) then
          f3(:nf,:4,i,j) = f3(:nf,:4,i,j) + cval*d3(:nf,:4,k,l)! (ij|lk)
          f3(:nf,:4,k,l) = f3(:nf,:4,k,l) + cval*d3(:nf,:4,i,j)! (kl|ji)
          f3(:nf,:4,i,j) = f3(:nf,:4,i,j) + cval*d3(:nf,:4,l,k)! (ij|kl)
          f3(:nf,:4,l,k) = f3(:nf,:4,l,k) + cval*d3(:nf,:4,i,j)! (lk|ji)
          f3(:nf,:4,j,i) = f3(:nf,:4,j,i) + cval*d3(:nf,:4,k,l)! (ji|lk)
          f3(:nf,:4,k,l) = f3(:nf,:4,k,l) + cval*d3(:nf,:4,j,i)! (kl|ij)
          f3(:nf,:4,j,i) = f3(:nf,:4,j,i) + cval*d3(:nf,:4,l,k)! (ji|kl)
          f3(:nf,:4,l,k) = f3(:nf,:4,l,k) + cval*d3(:nf,:4,j,i)! (lk|ij)

          f3(:nf,:7,i,k) = f3(:nf,:7,i,k) - xval*d3(:nf,:7,j,l)
          f3(:nf,:7,k,i) = f3(:nf,:7,k,i) - xval*d3(:nf,:7,l,j)
          f3(:nf,:7,i,l) = f3(:nf,:7,i,l) - xval*d3(:nf,:7,j,k)
          f3(:nf,:7,l,i) = f3(:nf,:7,l,i) - xval*d3(:nf,:7,k,j)
          f3(:nf,:7,j,k) = f3(:nf,:7,j,k) - xval*d3(:nf,:7,i,l)
          f3(:nf,:7,k,j) = f3(:nf,:7,k,j) - xval*d3(:nf,:7,l,i)
          f3(:nf,:7,j,l) = f3(:nf,:7,j,l) - xval*d3(:nf,:7,i,k)
          f3(:nf,:7,l,j) = f3(:nf,:7,l,j) - xval*d3(:nf,:7,k,i)
        else if (this%cur_pass==2) then
          f3(1:nf,7,i,k) = f3(1:nf,7,i,k) - xval*d3(1:nf,7,j,l)
          f3(1:nf,7,k,i) = f3(1:nf,7,k,i) - xval*d3(1:nf,7,l,j)
          f3(1:nf,7,i,l) = f3(1:nf,7,i,l) - xval*d3(1:nf,7,j,k)
          f3(1:nf,7,l,i) = f3(1:nf,7,l,i) - xval*d3(1:nf,7,k,j)
          f3(1:nf,7,j,k) = f3(1:nf,7,j,k) - xval*d3(1:nf,7,i,l)
          f3(1:nf,7,k,j) = f3(1:nf,7,k,j) - xval*d3(1:nf,7,l,i)
          f3(1:nf,7,j,l) = f3(1:nf,7,j,l) - xval*d3(1:nf,7,i,k)
          f3(1:nf,7,l,j) = f3(1:nf,7,l,j) - xval*d3(1:nf,7,k,i)
        end if

      end do
    end associate

    buf%ncur = 0

  end subroutine

subroutine int2_umrsf_data_t_update(this, buf)

  use io_constants, only: iw
    
  implicit none
  logical :: debug_mode
  class(int2_umrsf_data_t), intent(inout) :: this
  type(int2_storage_t), intent(inout) :: buf
  integer :: i, j, k, l, n
  real(kind=dp) :: val, xval, cval
  integer :: mythread, gpu_ierr

  mythread = buf%thread_id

  debug_mode = .True.

  if (.not.this%tamm_dancoff) return

  if (gpu_backend_metc_enabled()) then
    ! Resident hot path (UMRSF): borrow workspace pointers, per-thread scratch
    ! and f3 slice; no per-flush cudaMalloc/cudaFree.
    call gpu_backend_metc_try_contract(this%nfocks, ubound(this%d3,2), &
                                       ubound(this%d3,3), this%nthreads, &
                                       this%cur_pass, buf%buf_size, this%d3, &
                                       buf%thread_id - 1, buf%ids(:,1:buf%ncur), &
                                       buf%ints(1:buf%ncur), buf%ncur, &
                                       this%scale_exchange, this%scale_coulomb, &
                                       .true., gpu_ierr)
    if (gpu_ierr == 0) then
      buf%ncur = 0
      return
    end if
  end if

  associate ( f3 => this%f3(:,:,:,:,mythread), &
              d3 => this%d3, &
              nf => this%nfocks &
    )


    do n = 1, buf%ncur
      i = buf%ids(1,n)
      j = buf%ids(2,n)
      k = buf%ids(3,n)
      l = buf%ids(4,n)
      val = buf%ints(n)


      xval = val * this%scale_exchange
      cval = val * this%scale_coulomb

      if (this%cur_pass==1) then
        ! --- Coulomb-like updates (was :4 -> now :8 (alpha/beta pairs)) ---
        f3(:nf,1:8,i,j) = f3(:nf,1:8,i,j) + cval*d3(:nf,1:8,k,l)   ! (ij|lk)
        f3(:nf,1:8,k,l) = f3(:nf,1:8,k,l) + cval*d3(:nf,1:8,i,j)   ! (kl|ji)
        f3(:nf,1:8,i,j) = f3(:nf,1:8,i,j) + cval*d3(:nf,1:8,l,k)   ! (ij|kl)
        f3(:nf,1:8,l,k) = f3(:nf,1:8,l,k) + cval*d3(:nf,1:8,i,j)   ! (lk|ji)
        f3(:nf,1:8,j,i) = f3(:nf,1:8,j,i) + cval*d3(:nf,1:8,k,l)   ! (ji|lk)
        f3(:nf,1:8,k,l) = f3(:nf,1:8,k,l) + cval*d3(:nf,1:8,j,i)   ! (kl|ij)
        f3(:nf,1:8,j,i) = f3(:nf,1:8,j,i) + cval*d3(:nf,1:8,l,k)   ! (ji|kl)
        f3(:nf,1:8,l,k) = f3(:nf,1:8,l,k) + cval*d3(:nf,1:8,j,i)   ! (lk|ij)
        ! --- Exchange-like updates (was :7 -> now :11 (all types incl alpha/beta)) ---
        f3(:nf,1:8,i,k) = f3(:nf,1:8,i,k) - xval*d3(:nf,1:8,j,l) ! (ij|lk)
        f3(:nf,1:8,k,i) = f3(:nf,1:8,k,i) - xval*d3(:nf,1:8,l,j) ! (kl|ji)
        f3(:nf,1:8,i,l) = f3(:nf,1:8,i,l) - xval*d3(:nf,1:8,j,k) ! (ij|kl)
        f3(:nf,1:8,l,i) = f3(:nf,1:8,l,i) - xval*d3(:nf,1:8,k,j) ! (lk|ji)
        f3(:nf,1:8,j,k) = f3(:nf,1:8,j,k) - xval*d3(:nf,1:8,i,l) ! (ji|lk)
        f3(:nf,1:8,k,j) = f3(:nf,1:8,k,j) - xval*d3(:nf,1:8,l,i) ! (kl|ij)
        f3(:nf,1:8,j,l) = f3(:nf,1:8,j,l) - xval*d3(:nf,1:8,i,k) ! (ji|kl)
        f3(:nf,1:8,l,j) = f3(:nf,1:8,l,j) - xval*d3(:nf,1:8,k,i) ! (lk|ij)

        f3(:nf,9:10,i,l) = f3(:nf,9:10,i,l) - xval*d3(:nf,9:10,k,j) ! 
        f3(:nf,9:10,l,i) = f3(:nf,9:10,l,i) - xval*d3(:nf,9:10,j,k) ! 
        f3(:nf,9:10,k,j) = f3(:nf,9:10,k,j) - xval*d3(:nf,9:10,i,l) ! 
        f3(:nf,9:10,j,k) = f3(:nf,9:10,j,k) - xval*d3(:nf,9:10,l,i) ! 
        f3(:nf,9:10,i,k) = f3(:nf,9:10,i,k) - xval*d3(:nf,9:10,l,j) ! 
        f3(:nf,9:10,k,i) = f3(:nf,9:10,k,i) - xval*d3(:nf,9:10,j,l) ! 
        f3(:nf,9:10,l,j) = f3(:nf,9:10,l,j) - xval*d3(:nf,9:10,i,k) ! 
        f3(:nf,9:10,j,l) = f3(:nf,9:10,j,l) - xval*d3(:nf,9:10,k,i) ! 

        f3(1:nf,11,i,k) = f3(1:nf,11,i,k) - xval*d3(1:nf,11,j,l)
        f3(1:nf,11,k,i) = f3(1:nf,11,k,i) - xval*d3(1:nf,11,l,j)
        f3(1:nf,11,i,l) = f3(1:nf,11,i,l) - xval*d3(1:nf,11,j,k)
        f3(1:nf,11,l,i) = f3(1:nf,11,l,i) - xval*d3(1:nf,11,k,j)
        f3(1:nf,11,j,k) = f3(1:nf,11,j,k) - xval*d3(1:nf,11,i,l)
        f3(1:nf,11,k,j) = f3(1:nf,11,k,j) - xval*d3(1:nf,11,l,i)
        f3(1:nf,11,j,l) = f3(1:nf,11,j,l) - xval*d3(1:nf,11,i,k)
        f3(1:nf,11,l,j) = f3(1:nf,11,l,j) - xval*d3(1:nf,11,k,i)

      else if (this%cur_pass==2) then
        ! In pass 2 only agdlr was updated in scalar version (col 7).
        f3(1:nf,11,i,k) = f3(1:nf,11,i,k) - xval*d3(1:nf,11,j,l)
        f3(1:nf,11,k,i) = f3(1:nf,11,k,i) - xval*d3(1:nf,11,l,j)
        f3(1:nf,11,i,l) = f3(1:nf,11,i,l) - xval*d3(1:nf,11,j,k)
        f3(1:nf,11,l,i) = f3(1:nf,11,l,i) - xval*d3(1:nf,11,k,j)
        f3(1:nf,11,j,k) = f3(1:nf,11,j,k) - xval*d3(1:nf,11,i,l)
        f3(1:nf,11,k,j) = f3(1:nf,11,k,j) - xval*d3(1:nf,11,l,i)
        f3(1:nf,11,j,l) = f3(1:nf,11,j,l) - xval*d3(1:nf,11,i,k)
        f3(1:nf,11,l,j) = f3(1:nf,11,l,j) - xval*d3(1:nf,11,k,i)
        ! Here agdlr is column 11 (spin-independent), update only that.
      end if

    end do

!  if (debug_mode ) then
!    write(iw,*) 'UPDATE'
!    write(iw,*) 'adco2a 1 1-5 1', f3(1,7,1:5,1)
!    write(iw,*) 'agdlr 1 1 1-5', f3(1,11,1,1:5) 
!    write(iw,*) this%scale_exchange, this%scale_coulomb !xval, cval
!  endif

  end associate

  buf%ncur = 0

end subroutine

