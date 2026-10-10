!===============================================================================
!> @brief Route-C response bridge (APIV2.md R3c, libintRot repo).
!>
!> Optional external supplier of the TDDFT/SF/MRSF response J/K builds.
!> If $OQP_ROUTEC_LIB names a dylib exposing the APIV2 symbols
!>
!>   int  routec_jkm_init(const int* nbf);                    // 0 = ok
!>   void routec_jkm_apply(const double* x, const int* nbf, const int* nvec,
!>                         double* jm, double* km, int* info); // info 0 = ok
!>
!> then int2_run diverts the int2_td_data_t / int2_tdgrd_data_t /
!> int2_mrsf_data_t consumers (EXACT dynamic types only —
!> int2_umrsf_data_t / int2_rpagrd_data_t fall through to the native
!> quartet loop) to batched RAW J(M)/K(M) contractions
!> of the trial densities, with
!>
!>   J(M)_ab = (ab|kl) M_kl ,  K(M)_ab = (ak|bl) M_kl   (M non-symmetric).
!>
!> All consumer prefactors are applied HERE, from the verified Stage-0
!> tables (OPENQP_BRIDGE.md "Stage 0 RESULTS"), in the PRE-parallel_stop,
!> thread-slot-1 deposit convention:
!>   TD RPA  : apb deposit = half of [2c*J(Ms) - x*K(Ms)], Ms = d2 + d2^T
!>             (parallel_stop applies A := A + A^T, diagonal included);
!>             amb deposit = x*[K(d2)^T - K(d2)]            (full value)
!>   TD TDA  : amb deposit = 2c*J(d2)*[tdc] - x*K(d2)       (full value)
!>   TDGRD   : (d2 = alpha/beta pair, G-5 verified) apb_s deposit = half of
!>             [c*J(Ms_1+Ms_2) - x*K(Ms_s)]; amb slot 1 only, as TD RPA
!>   MRSF    : f3 slot m deposit = c*J(D_m)*[m<=4] - x*K(D_m) (full value)
!> using the single-apply identities J(Ms) = 2 J(d2), K(Ms) = K(d2)+K(d2)^T.
!>
!> Compute on rank 0 only (other ranks deposit nothing; parallel_stop
!> allreduces). CAM passes never reach this hook: int2_run_cam calls
!> run_generic directly, so attenuated/short-range builds stay native.
!> Entirely inert when $OQP_ROUTEC_LIB is unset or lacks the symbols; falls
!> back natively whenever init/apply report failure.
!>
!> APIV2.1 (low-rank fast path): when the dylib ALSO exports
!>
!>   void routec_jkm_apply_lr(const double* u, const double* v,
!>        const int* nbf, const int* nout, const int* ranks,
!>        double* jm, double* km, int* info);   // M_o = sum_r u_r v_r^T
!>
!> and the MRSF consumer carries the mrsfcbc factor table (d3fac), slots
!> 1-6 are routed through the low-rank apply (slots 1-4 rank-1, 5-6
!> rank-2) and only slot 7 stays dense — same deposits, ~nbf/2x fewer
!> FLOPs on 6 of 7 slots. OQP_ROUTEC_LR=0 forces the dense branch (the
!> validation referee); any lr failure also falls back to it.
!>
!> This is a SUBMODULE of int2_compute so it may use tdhf_lib /
!> tdhf_mrsf_lib (which themselves use int2_compute) without a module cycle.
!===============================================================================
submodule (int2_compute) int2_routec
  use, intrinsic :: iso_c_binding
  implicit none

  abstract interface
    function routec_jkm_init_i(nbf) result(ierr) bind(C)
      import :: c_int
      integer(c_int), intent(in) :: nbf
      integer(c_int) :: ierr
    end function
    subroutine routec_jkm_apply_i(x, nbf, nvec, jm, km, info) bind(C)
      import :: c_ptr, c_int
      type(c_ptr), value :: x, jm, km
      integer(c_int), intent(in) :: nbf, nvec
      integer(c_int), intent(out) :: info
    end subroutine
    subroutine routec_jkm_apply_lr_i(u, v, nbf, nout, ranks, jm, km, info) &
        bind(C)
      import :: c_ptr, c_int
      type(c_ptr), value :: u, v, jm, km
      integer(c_int), intent(in) :: nbf, nout
      integer(c_int), intent(in) :: ranks(*)
      integer(c_int), intent(out) :: info
    end subroutine
  end interface

  interface
    function rtc_dlopen(file, mode) bind(C, name="dlopen")
      import :: c_ptr, c_char, c_int
      character(kind=c_char), intent(in) :: file(*)
      integer(c_int), value :: mode
      type(c_ptr) :: rtc_dlopen
    end function
    function rtc_dlsym(handle, name) bind(C, name="dlsym")
      import :: c_ptr, c_char, c_funptr
      type(c_ptr), value :: handle
      character(kind=c_char), intent(in) :: name(*)
      type(c_funptr) :: rtc_dlsym
    end function
  end interface

  integer(c_int), parameter :: RTC_RTLD_NOW = 2_c_int
  logical :: rtc_tried = .false.
  integer :: rtc_nbf_ok = -1
  logical :: rtc_lr_on = .true.        ! OQP_ROUTEC_LR=0 forces dense branch
  procedure(routec_jkm_init_i), pointer :: rtc_init => null()
  procedure(routec_jkm_apply_i), pointer :: rtc_apply => null()
  procedure(routec_jkm_apply_lr_i), pointer :: rtc_apply_lr => null()

contains

  !> One-time dlopen/dlsym of $OQP_ROUTEC_LIB, then per-nbf session init.
  logical function rtc_session(nbf)
    integer, intent(in) :: nbf
    type(c_ptr) :: h
    type(c_funptr) :: fp
    character(len=4096) :: path
    integer :: plen, stat
    integer(c_int) :: nbf_c, ierr
    rtc_session = .false.
    if (.not. rtc_tried) then
      rtc_tried = .true.
      call get_environment_variable("OQP_ROUTEC_LIB", path, plen, stat)
      if (stat /= 0 .or. plen == 0) return
      h = rtc_dlopen(trim(path)//c_null_char, RTC_RTLD_NOW)
      if (.not. c_associated(h)) return
      fp = rtc_dlsym(h, "routec_jkm_init"//c_null_char)
      if (.not. c_associated(fp)) return
      call c_f_procpointer(fp, rtc_init)
      fp = rtc_dlsym(h, "routec_jkm_apply"//c_null_char)
      if (.not. c_associated(fp)) then
        rtc_init => null()
        return
      end if
      call c_f_procpointer(fp, rtc_apply)
      ! v2.1 low-rank apply is OPTIONAL: absent => dense applies only
      fp = rtc_dlsym(h, "routec_jkm_apply_lr"//c_null_char)
      if (c_associated(fp)) call c_f_procpointer(fp, rtc_apply_lr)
      call get_environment_variable("OQP_ROUTEC_LR", path, plen, stat)
      if (stat == 0 .and. plen > 0) rtc_lr_on = (path(1:1) /= '0')
    end if
    if (.not. associated(rtc_apply)) return
    if (rtc_nbf_ok /= nbf) then
      nbf_c = int(nbf, c_int)
      ierr = rtc_init(nbf_c)
      if (ierr /= 0_c_int) return
      rtc_nbf_ok = nbf
    end if
    rtc_session = .true.
  end function rtc_session

  !> Hook called from int2_run (non-CAM path only). .true. = handled.
  module function int2_routec_response(this, int2_consumer) result(done)
    use tdhf_lib, only: int2_td_data_t, int2_tdgrd_data_t
    use tdhf_mrsf_lib, only: int2_mrsf_data_t
    class(int2_compute_t), intent(inout) :: this
    class(int2_compute_data_t), intent(inout) :: int2_consumer
    logical :: done
    done = .false.
    if (this%attenuated) return
    select type (c => int2_consumer)
    type is (int2_tdgrd_data_t)
      done = rtc_tdgrd(this, c)
    type is (int2_td_data_t)
      done = rtc_td(this, c)
    type is (int2_mrsf_data_t)
      done = rtc_mrsf(this, c)
    end select
  end function int2_routec_response

  !> int2_tdgrd_data_t (Z-vector / gradient response builds): d2(nbf,nbf,2)
  !> alpha/beta trial pair. Verified table (G-5 phase A(i), probe
  !> g5_probeA_zvec.py, relaxed density + W rebuilt to ~2e-8 with
  !> wrong-prefactor controls failing at 1e-2..1e0):
  !>   apb_s = c*J(Ms_1 + Ms_2) - x*K(Ms_s),  Ms_s := d2_s + d2_s^T
  !> (Coulomb of the spin-SUMMED symmetrized density deposited in BOTH
  !> slots; exchange spin-diagonal). Deposit halves, A := A + A^T applied
  !> in parallel_stop. amb (int_amb, written for slot 1 only in the native
  !> update) = x*[K(d2_1)^T - K(d2_1)], full value, not symmetrized.
  logical function rtc_tdgrd(driver, dat)
    use tdhf_lib, only: int2_tdgrd_data_t
    class(int2_compute_t), intent(inout) :: driver
    type(int2_tdgrd_data_t), intent(inout) :: dat
    integer :: nbf, s
    integer(c_int) :: nbf_c, nvec_c, info
    real(c_double), allocatable, target :: x(:,:,:), jm(:,:,:), km(:,:,:)
    real(kind=dp) :: cc, xx
    rtc_tdgrd = .false.
    if (size(dat%d2, 3) /= 2) return     ! native update hardcodes 2 slots
    if (dat%tamm_dancoff) return         ! no TDA branch in the native update
    nbf = driver%basis%nbf
    if (.not. rtc_session(nbf)) return
    call dat%parallel_start(driver%basis, 1)
    if (driver%pe%rank == 0) then
      allocate(x(nbf,nbf,2), jm(nbf,nbf,2), km(nbf,nbf,2))
      x = dat%d2(:,:,1:2)
      nbf_c = int(nbf, c_int)
      nvec_c = 2_c_int
      info = 1_c_int
      call rtc_apply(c_loc(x), nbf_c, nvec_c, c_loc(jm), c_loc(km), info)
      if (info /= 0_c_int) return   ! fall back natively
      cc = dat%scale_coulomb
      xx = dat%scale_exchange
      if (dat%int_apb) then
        do s = 1, 2   ! half of c*J(Mtot) - x*K(Ms_s); A+A^T downstream
          dat%apb(:,:,s,1) = cc*(jm(:,:,1) + jm(:,:,2)) &
                           - 0.5_dp*xx*(km(:,:,s) + transpose(km(:,:,s)))
        end do
      end if
      if (dat%int_amb) &
        dat%amb(:,:,1,1) = xx*(transpose(km(:,:,1)) - km(:,:,1))
    end if
    call dat%pe%init(driver%pe%comm, driver%pe%use_mpi)
    call dat%parallel_stop()
    rtc_tdgrd = .true.
    call rtc_note_once(1, 'tdgrd')
  end function rtc_tdgrd

  !> One-time stderr notice per consumer kind (positive evidence which
  !> seams fired in a given run).
  subroutine rtc_note_once(kind, name)
    use, intrinsic :: iso_fortran_env, only: error_unit
    integer, intent(in) :: kind
    character(len=*), intent(in) :: name
    logical, save :: noted(4) = .false.
    if (noted(kind)) return
    noted(kind) = .true.
    write(error_unit, '(a)') &
      ' [int2_routec] '//name//' consumer diverted to external J/K'
  end subroutine rtc_note_once

  !> int2_td_data_t: d2(nbf,nbf,nvec) trial densities, one J+K apply of the
  !> raw d2 batch; deposits per the Stage-0 table (see header).
  logical function rtc_td(driver, dat)
    use tdhf_lib, only: int2_td_data_t
    class(int2_compute_t), intent(inout) :: driver
    type(int2_td_data_t), intent(inout) :: dat
    integer :: nbf, nvec, v
    integer(c_int) :: nbf_c, nvec_c, info
    real(c_double), allocatable, target :: x(:,:,:), jm(:,:,:), km(:,:,:)
    real(kind=dp) :: cc, xx
    rtc_td = .false.
    nbf = driver%basis%nbf
    nvec = size(dat%d2, 3)
    if (nvec < 1) return
    if (.not. rtc_session(nbf)) return
    call dat%parallel_start(driver%basis, 1)
    if (driver%pe%rank == 0) then
      allocate(x(nbf,nbf,nvec), jm(nbf,nbf,nvec), km(nbf,nbf,nvec))
      x = dat%d2(:,:,1:nvec)
      nbf_c = int(nbf, c_int)
      nvec_c = int(nvec, c_int)
      info = 1_c_int
      call rtc_apply(c_loc(x), nbf_c, nvec_c, c_loc(jm), c_loc(km), info)
      if (info /= 0_c_int) return   ! fall back natively (arrays re-zeroed there)
      cc = dat%scale_coulomb
      xx = dat%scale_exchange
      do v = 1, nvec
        if (dat%tamm_dancoff) then
          dat%amb(:,:,v,1) = -xx*km(:,:,v)
          if (dat%tamm_dancoff_coulomb) &
            dat%amb(:,:,v,1) = dat%amb(:,:,v,1) + 2.0_dp*cc*jm(:,:,v)
        else
          if (dat%int_apb) &   ! half of 2c*J(Ms) - x*K(Ms); A+A^T downstream
            dat%apb(:,:,v,1) = 2.0_dp*cc*jm(:,:,v) &
                             - 0.5_dp*xx*(km(:,:,v) + transpose(km(:,:,v)))
          if (dat%int_amb) &
            dat%amb(:,:,v,1) = xx*(transpose(km(:,:,v)) - km(:,:,v))
        end if
      end do
    end if
    call dat%pe%init(driver%pe%comm, driver%pe%use_mpi)
    call dat%parallel_stop()
    rtc_td = .true.
    call rtc_note_once(2, 'td')
  end function rtc_td

  !> int2_mrsf_data_t: d3(nvec,7,nbf,nbf), vector index FIRST. Preferred
  !> branch (APIV2.1): slots 1-6 through the low-rank apply from the
  !> mrsfcbc factor table dat%d3fac, slot 7 dense (rtc_mrsf_lr). Fallback
  !> branch (also the validation referee, forced by OQP_ROUTEC_LR=0):
  !> gather all 7 slots to X(nbf,nbf,nvec*7), one dense batched apply,
  !> scatter back to f3(nvec,7,:,:,1). Both deposit identical RAW values.
  logical function rtc_mrsf(driver, dat)
    use tdhf_mrsf_lib, only: int2_mrsf_data_t
    class(int2_compute_t), intent(inout) :: driver
    type(int2_mrsf_data_t), intent(inout) :: dat
    integer :: nbf, nvec, nmat, v, m, ib
    integer(c_int) :: nbf_c, nvec_c, info
    real(c_double), allocatable, target :: x(:,:,:), jm(:,:,:), km(:,:,:)
    real(kind=dp) :: cc, xx
    logical :: lr_done
    rtc_mrsf = .false.
    if (.not. dat%tamm_dancoff) return   ! native update is a no-op guard too
    nbf = driver%basis%nbf
    nvec = size(dat%d3, 1)
    nmat = size(dat%d3, 2)
    if (nvec < 1 .or. nmat /= 7) return
    if (.not. rtc_session(nbf)) return
    call dat%parallel_start(driver%basis, 1)
    if (driver%pe%rank == 0) then
      lr_done = .false.
      if (rtc_lr_on .and. associated(rtc_apply_lr) &
          .and. associated(dat%d3fac)) then
        if (size(dat%d3fac,1) == nbf .and. size(dat%d3fac,2) == 8 &
            .and. size(dat%d3fac,3) == nvec) lr_done = rtc_mrsf_lr(dat, nbf, nvec)
      end if
      if (.not. lr_done) then
        allocate(x(nbf,nbf,nvec*nmat), jm(nbf,nbf,nvec*nmat), km(nbf,nbf,nvec*nmat))
        do m = 1, nmat
          do v = 1, nvec
            ib = (m-1)*nvec + v
            x(:,:,ib) = dat%d3(v,m,:,:)
          end do
        end do
        nbf_c = int(nbf, c_int)
        nvec_c = int(nvec*nmat, c_int)
        info = 1_c_int
        call rtc_apply(c_loc(x), nbf_c, nvec_c, c_loc(jm), c_loc(km), info)
        if (info /= 0_c_int) return
        cc = dat%scale_coulomb
        xx = dat%scale_exchange
        do m = 1, nmat
          do v = 1, nvec
            ib = (m-1)*nvec + v
            if (m <= 4) then   ! Coulomb loop covers slots 1:4 only
              dat%f3(v,m,:,:,1) = cc*jm(:,:,ib) - xx*km(:,:,ib)
            else               ! slots 5-7 exchange-only
              dat%f3(v,m,:,:,1) = -xx*km(:,:,ib)
            end if
          end do
        end do
      end if
    end if
    call dat%pe%init(driver%pe%comm, driver%pe%use_mpi)
    call dat%parallel_stop()
    rtc_mrsf = .true.
    call rtc_note_once(3, 'mrsf')
  end function rtc_mrsf

  !> APIV2.1 low-rank MRSF branch (rank 0 only; caller did parallel_start).
  !> d3fac(nbf,8,v) columns are u1,v1,u2,v2,u3,v3,u4,v4 from mrsfcbc, with
  !> slot densities M_m = u_m v_m^T (m=1..4) and the rank-2 combinations
  !>   M5 = v1 u2^T - v2 u1^T ,   M6 = v4 u3^T - v3 u4^T
  !> (minus signs folded into the u-side factor columns below). Three
  !> applies, deposits identical to the dense branch (Stage-0 table):
  !>   A) slots 1-4: 4*nvec rank-1 outputs, jm+km, f3 = c*J - x*K
  !>   B) slots 5-6: 2*nvec rank-2 outputs, km only, f3 = -x*K
  !>   C) slot 7   : dense apply of d3(:,7,:,:), km only, f3 = -x*K
  !> Returns .false. on any engine failure; the caller then redoes the
  !> work through the dense branch (engine state is untouched on failure).
  logical function rtc_mrsf_lr(dat, nbf, nvec)
    use tdhf_mrsf_lib, only: int2_mrsf_data_t
    type(int2_mrsf_data_t), intent(inout) :: dat
    integer, intent(in) :: nbf, nvec
    integer :: v, m, o
    integer(c_int) :: nbf_c, n_c, info
    integer(c_int), allocatable :: ranks(:)
    real(c_double), allocatable, target :: uu(:,:), vv(:,:)
    real(c_double), allocatable, target :: jm(:,:,:), km(:,:,:), x7(:,:,:)
    real(kind=dp) :: cc, xx
    rtc_mrsf_lr = .false.
    cc = dat%scale_coulomb
    xx = dat%scale_exchange
    nbf_c = int(nbf, c_int)

    ! ---- A) slots 1-4, rank-1: output o = (m-1)*nvec + v
    allocate(uu(nbf,4*nvec), vv(nbf,4*nvec), ranks(4*nvec), &
             jm(nbf,nbf,4*nvec), km(nbf,nbf,4*nvec))
    ranks = 1_c_int
    do m = 1, 4
      do v = 1, nvec
        o = (m-1)*nvec + v
        uu(:,o) = dat%d3fac(:,2*m-1,v)
        vv(:,o) = dat%d3fac(:,2*m,  v)
      end do
    end do
    n_c = int(4*nvec, c_int)
    info = 1_c_int
    call rtc_apply_lr(c_loc(uu), c_loc(vv), nbf_c, n_c, ranks, &
                      c_loc(jm), c_loc(km), info)
    if (info /= 0_c_int) return
    do m = 1, 4
      do v = 1, nvec
        o = (m-1)*nvec + v
        dat%f3(v,m,:,:,1) = cc*jm(:,:,o) - xx*km(:,:,o)
      end do
    end do
    deallocate(uu, vv, ranks, jm, km)

    ! ---- B) slots 5-6, rank-2: outputs v (m=5) and nvec+v (m=6),
    !      factor columns 2*o-1, 2*o per output
    allocate(uu(nbf,4*nvec), vv(nbf,4*nvec), ranks(2*nvec), &
             km(nbf,nbf,2*nvec))
    ranks = 2_c_int
    do v = 1, nvec
      o = v                      ! M5 = v1 u2^T - v2 u1^T
      uu(:,2*o-1) =  dat%d3fac(:,2,v)
      vv(:,2*o-1) =  dat%d3fac(:,3,v)
      uu(:,2*o)   = -dat%d3fac(:,4,v)
      vv(:,2*o)   =  dat%d3fac(:,1,v)
      o = nvec + v               ! M6 = v4 u3^T - v3 u4^T
      uu(:,2*o-1) =  dat%d3fac(:,8,v)
      vv(:,2*o-1) =  dat%d3fac(:,5,v)
      uu(:,2*o)   = -dat%d3fac(:,6,v)
      vv(:,2*o)   =  dat%d3fac(:,7,v)
    end do
    n_c = int(2*nvec, c_int)
    info = 1_c_int
    call rtc_apply_lr(c_loc(uu), c_loc(vv), nbf_c, n_c, ranks, &
                      c_null_ptr, c_loc(km), info)
    if (info /= 0_c_int) return
    do v = 1, nvec
      dat%f3(v,5,:,:,1) = -xx*km(:,:,v)
      dat%f3(v,6,:,:,1) = -xx*km(:,:,nvec+v)
    end do
    deallocate(uu, vv, ranks, km)

    ! ---- C) slot 7 dense, exchange-only
    allocate(x7(nbf,nbf,nvec), km(nbf,nbf,nvec))
    do v = 1, nvec
      x7(:,:,v) = dat%d3(v,7,:,:)
    end do
    n_c = int(nvec, c_int)
    info = 1_c_int
    call rtc_apply(c_loc(x7), nbf_c, n_c, c_null_ptr, c_loc(km), info)
    if (info /= 0_c_int) return
    do v = 1, nvec
      dat%f3(v,7,:,:,1) = -xx*km(:,:,v)
    end do
    rtc_mrsf_lr = .true.
    call rtc_note_once(4, 'mrsf low-rank (slots 1-6 factored, 7 dense)')
  end function rtc_mrsf_lr

end submodule int2_routec
