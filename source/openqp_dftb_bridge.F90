module openqp_dftb_bridge_mod
  use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int32_t, c_int64_t, c_null_char
  use c_interop, only: oqp_handle_t, oqp_handle_get_info
  use types, only: information
#ifdef OQP_ENABLE_OPENQP_DFTB
  use openqp_dftb_api, only: &
      openqp_dftb_run_ground_state, &
      openqp_dftb_run_ground_state_gradient, &
      openqp_dftb_run_mrsf_state_gradient, &
      openqp_dftb_run_mrsf_tddftb_response, &
      openqp_dftb_run_spin_flip_state_gradient, &
      openqp_dftb_run_spin_flip_tda_response, &
      openqp_dftb_run_tda_response
  use openqp_dftb_kinds, only: dftb_i8, dftb_dp
  use openqp_dftb_parameters, only: &
      openqp_dftb_read_parameter_file, &
      openqp_dftb_read_skf_directory
  use openqp_dftb_types, only: &
      openqp_dftb_atom_t, &
      openqp_dftb_options_t, &
      openqp_dftb_parameter_set_t, &
      openqp_dftb_result_t, &
      openqp_dftb_status_not_implemented, &
      openqp_dftb_status_ok
#endif
  implicit none
  private

  public :: oqp_dftb_state_gradient_C

contains

  subroutine oqp_dftb_state_gradient_C(c_handle, parameter_path, parameter_len, method_name, method_len, &
      root, n_roots, need_gradient, molecular_charge, scc_mixer, scc_history, max_scc_iterations, &
      response_max_iterations, response_max_subspace, scc_tolerance, scc_mixing, &
      scc_max_step, response_tolerance, spc, omega, cam_alpha, cam_beta, &
      reference_energy, state_energy, excitation_energy, spin_square, gradient, &
      status_message, status_message_len, status) &
      bind(C, name="oqp_dftb_state_gradient")
    type(oqp_handle_t) :: c_handle
    character(kind=c_char), intent(in) :: parameter_path(*)
    integer(c_int32_t), value :: parameter_len
    character(kind=c_char), intent(in) :: method_name(*)
    integer(c_int32_t), value :: method_len
    integer(c_int64_t), value :: root
    integer(c_int64_t), value :: n_roots
    integer(c_int64_t), value :: need_gradient
    integer(c_int64_t), value :: molecular_charge
    integer(c_int64_t), value :: scc_mixer
    integer(c_int64_t), value :: scc_history
    integer(c_int64_t), value :: max_scc_iterations
    integer(c_int64_t), value :: response_max_iterations
    integer(c_int64_t), value :: response_max_subspace
    real(c_double), value :: scc_tolerance
    real(c_double), value :: scc_mixing
    real(c_double), value :: scc_max_step
    real(c_double), value :: response_tolerance
    real(c_double), value :: spc
    real(c_double), value :: omega
    real(c_double), value :: cam_alpha
    real(c_double), value :: cam_beta
    real(c_double), intent(out) :: reference_energy
    real(c_double), intent(out) :: state_energy
    real(c_double), intent(out) :: excitation_energy
    real(c_double), intent(out) :: spin_square
    real(c_double), intent(out) :: gradient(3, *)
    character(kind=c_char), intent(out) :: status_message(*)
    integer(c_int32_t), value :: status_message_len
    integer(c_int64_t), intent(out) :: status

    reference_energy = 0.0_c_double
    state_energy = 0.0_c_double
    excitation_energy = 0.0_c_double
    spin_square = 0.0_c_double
    call write_c_string("", status_message, status_message_len)

#ifdef OQP_ENABLE_OPENQP_DFTB
    call run_openqp_dftb(c_handle, parameter_path, parameter_len, method_name, method_len, &
        root, n_roots, need_gradient, molecular_charge, scc_mixer, scc_history, max_scc_iterations, &
        response_max_iterations, response_max_subspace, scc_tolerance, scc_mixing, &
        scc_max_step, response_tolerance, spc, omega, cam_alpha, cam_beta, &
        reference_energy, state_energy, excitation_energy, spin_square, gradient, &
        status_message, status_message_len, status)
#else
    status = -9001_c_int64_t
    call write_c_string("OpenQP was built without ENABLE_OPENQP_DFTB", &
        status_message, status_message_len)
#endif
  end subroutine oqp_dftb_state_gradient_C

  subroutine write_c_string(text, chars, n_chars)
    character(len=*), intent(in) :: text
    character(kind=c_char), intent(out) :: chars(*)
    integer(c_int32_t), value :: n_chars
    integer :: i
    integer :: n

    if (n_chars <= 0) return
    do i = 1, int(n_chars)
      chars(i) = c_null_char
    enddo
    n = min(len_trim(text), max(0, int(n_chars) - 1))
    do i = 1, n
      chars(i) = text(i:i)
    enddo
  end subroutine write_c_string

#ifdef OQP_ENABLE_OPENQP_DFTB
  subroutine run_openqp_dftb(c_handle, parameter_path, parameter_len, method_name, method_len, &
      root, n_roots, need_gradient, molecular_charge, scc_mixer, scc_history, max_scc_iterations, &
      response_max_iterations, response_max_subspace, scc_tolerance, scc_mixing, &
      scc_max_step, response_tolerance, spc, omega, cam_alpha, cam_beta, &
      reference_energy, state_energy, excitation_energy, spin_square, gradient, &
      status_message, status_message_len, status)
    type(oqp_handle_t), intent(in) :: c_handle
    character(kind=c_char), intent(in) :: parameter_path(*)
    integer(c_int32_t), value :: parameter_len
    character(kind=c_char), intent(in) :: method_name(*)
    integer(c_int32_t), value :: method_len
    integer(c_int64_t), value :: root
    integer(c_int64_t), value :: n_roots
    integer(c_int64_t), value :: need_gradient
    integer(c_int64_t), value :: molecular_charge
    integer(c_int64_t), value :: scc_mixer
    integer(c_int64_t), value :: scc_history
    integer(c_int64_t), value :: max_scc_iterations
    integer(c_int64_t), value :: response_max_iterations
    integer(c_int64_t), value :: response_max_subspace
    real(c_double), value :: scc_tolerance
    real(c_double), value :: scc_mixing
    real(c_double), value :: scc_max_step
    real(c_double), value :: response_tolerance
    real(c_double), value :: spc
    real(c_double), value :: omega
    real(c_double), value :: cam_alpha
    real(c_double), value :: cam_beta
    real(c_double), intent(out) :: reference_energy
    real(c_double), intent(out) :: state_energy
    real(c_double), intent(out) :: excitation_energy
    real(c_double), intent(out) :: spin_square
    real(c_double), intent(out) :: gradient(3, *)
    character(kind=c_char), intent(out) :: status_message(*)
    integer(c_int32_t), value :: status_message_len
    integer(c_int64_t), intent(out) :: status

    type(information), pointer :: inf
    type(openqp_dftb_atom_t), allocatable :: atoms(:)
    type(openqp_dftb_parameter_set_t) :: parameters
    type(openqp_dftb_options_t) :: reference_options
    type(openqp_dftb_options_t) :: response_options
    type(openqp_dftb_result_t) :: result
    integer(dftb_i8), allocatable :: unique_atomic_numbers(:)
    real(dftb_dp), allocatable :: dftb_gradient(:,:)
    character(len=:), allocatable :: parameter
    character(len=:), allocatable :: method
    character(len=:), allocatable :: message
    integer(dftb_i8) :: dftb_status
    integer :: iatom
    integer :: natom
    logical :: compute_gradient

    inf => oqp_handle_get_info(c_handle)
    if (.not. allocated(inf%atoms%xyz) .or. .not. allocated(inf%atoms%zn)) then
      status = -9002_c_int64_t
      call write_c_string("OpenQP molecule handle is missing atoms or coordinates", &
          status_message, status_message_len)
      return
    endif

    natom = size(inf%atoms%zn)
    gradient(:, 1:natom) = 0.0_c_double
    allocate(atoms(natom))
    do iatom = 1, natom
      atoms(iatom)%atomic_number = int(nint(inf%atoms%zn(iatom)), dftb_i8)
      atoms(iatom)%xyz_bohr = real(inf%atoms%xyz(:, iatom), dftb_dp)
    enddo

    parameter = c_string(parameter_path, parameter_len)
    method = lowercase(c_string(method_name, method_len))
    compute_gradient = need_gradient /= 0_c_int64_t

    if (index(trim(parameter), ".opdftb") > 0) then
      call openqp_dftb_read_parameter_file(trim(parameter), parameters, dftb_status, message)
    else
      call collect_unique_atomic_numbers(atoms, unique_atomic_numbers)
      call openqp_dftb_read_skf_directory(trim(parameter), unique_atomic_numbers, parameters, dftb_status, message)
    endif
    if (dftb_status /= openqp_dftb_status_ok) then
      status = int(dftb_status, c_int64_t)
      call write_c_string(message, status_message, status_message_len)
      return
    endif

    reference_options = openqp_dftb_options_t()
    reference_options%scc = .true.
    reference_options%molecular_charge = int(molecular_charge, dftb_i8)
    reference_options%scc_mixer = int(scc_mixer, dftb_i8)
    reference_options%scc_mixer_history = int(scc_history, dftb_i8)
    reference_options%max_scc_iterations = int(max_scc_iterations, dftb_i8)
    reference_options%response_max_iterations = int(response_max_iterations, dftb_i8)
    reference_options%response_max_subspace = int(response_max_subspace, dftb_i8)
    reference_options%n_roots = max(int(n_roots, dftb_i8), max(1_dftb_i8, int(root, dftb_i8)))
    reference_options%scc_tolerance = real(scc_tolerance, dftb_dp)
    reference_options%scc_mixing = real(scc_mixing, dftb_dp)
    reference_options%scc_mixer_max_step = real(scc_max_step, dftb_dp)
    reference_options%response_tolerance = real(response_tolerance, dftb_dp)
    reference_options%spc_coco = real(spc, dftb_dp)
    reference_options%spc_ovov = real(spc, dftb_dp)
    reference_options%spc_coov = real(spc, dftb_dp)
    reference_options%range_separation_omega = real(omega, dftb_dp)
    reference_options%cam_alpha = real(cam_alpha, dftb_dp)
    reference_options%cam_beta = real(cam_beta, dftb_dp)

    select case (trim(method))
    case ("ground", "dftb")
      reference_options%spin_polarized = .false.
      reference_options%spin_complete = .false.
      reference_options%mrsf = .false.
      reference_options%multiplicity = 1_dftb_i8
      call openqp_dftb_run_ground_state(atoms, reference_options, parameters, result, dftb_status, message)
      if (compute_gradient .and. dftb_status == openqp_dftb_status_ok) then
        call openqp_dftb_run_ground_state_gradient(atoms, reference_options, parameters, result, &
            dftb_gradient, dftb_status, message)
      endif
      state_energy = result%total_energy
    case ("ground_noscc", "noscc", "dftb0")
      reference_options%scc = .false.
      reference_options%spin_polarized = .false.
      reference_options%spin_complete = .false.
      reference_options%mrsf = .false.
      reference_options%multiplicity = 1_dftb_i8
      call openqp_dftb_run_ground_state(atoms, reference_options, parameters, result, dftb_status, message)
      if (compute_gradient .and. dftb_status == openqp_dftb_status_ok) then
        call openqp_dftb_run_ground_state_gradient(atoms, reference_options, parameters, result, &
            dftb_gradient, dftb_status, message)
      endif
      state_energy = result%total_energy
    case ("tddftb", "tda")
      reference_options%spin_polarized = .false.
      reference_options%spin_complete = .false.
      reference_options%mrsf = .false.
      reference_options%multiplicity = 1_dftb_i8
      call openqp_dftb_run_ground_state(atoms, reference_options, parameters, result, dftb_status, message)
      if (dftb_status == openqp_dftb_status_ok) then
        response_options = reference_options
        call openqp_dftb_run_tda_response(response_options, result, dftb_status, message)
      endif
      if (dftb_status == openqp_dftb_status_ok) then
        call check_response_root(root, result%n_roots, dftb_status)
      endif
      if (dftb_status == openqp_dftb_status_ok) then
        state_energy = result%total_energy + result%excitation_energy(root)
        excitation_energy = result%excitation_energy(root)
        spin_square = result%spin_square(root)
        if (compute_gradient) then
          dftb_status = openqp_dftb_status_not_implemented
          message = "OpenQP-DFTB TDDFTB/TDA state gradient is not exported by the current API"
        endif
      endif
    case ("sf", "sftddftb", "sf-tddftb")
      call set_spin_flip_reference(reference_options)
      call openqp_dftb_run_ground_state(atoms, reference_options, parameters, result, dftb_status, message)
      if (dftb_status == openqp_dftb_status_ok) then
        response_options = reference_options
        response_options%multiplicity = 1_dftb_i8
        call openqp_dftb_run_spin_flip_tda_response(response_options, result, dftb_status, message)
      endif
      if (dftb_status == openqp_dftb_status_ok) then
        call check_response_root(root, result%n_roots, dftb_status)
      endif
      if (dftb_status == openqp_dftb_status_ok) then
        state_energy = result%total_energy + result%excitation_energy(root)
        excitation_energy = result%excitation_energy(root)
        spin_square = result%spin_square(root)
        if (compute_gradient) then
          call openqp_dftb_run_spin_flip_state_gradient(atoms, response_options, parameters, result, &
              int(root), dftb_gradient, dftb_status, message)
        endif
      endif
    case ("mrsf", "mrsftddftb", "mrsf-tddftb")
      call set_spin_flip_reference(reference_options)
      call openqp_dftb_run_ground_state(atoms, reference_options, parameters, result, dftb_status, message)
      if (dftb_status == openqp_dftb_status_ok) then
        response_options = reference_options
        response_options%mrsf = .true.
        response_options%spin_complete = .true.
        response_options%multiplicity = 1_dftb_i8
        call openqp_dftb_run_mrsf_tddftb_response(response_options, result, dftb_status, message)
      endif
      if (dftb_status == openqp_dftb_status_ok) then
        call check_response_root(root, result%n_roots, dftb_status)
      endif
      if (dftb_status == openqp_dftb_status_ok) then
        state_energy = result%total_energy + result%excitation_energy(root)
        excitation_energy = result%excitation_energy(root)
        spin_square = result%spin_square(root)
        if (compute_gradient) then
          call openqp_dftb_run_mrsf_state_gradient(atoms, response_options, parameters, result, &
              int(root), dftb_gradient, dftb_status, message)
        endif
      endif
    case default
      dftb_status = -9003_dftb_i8
      message = "Unknown OpenQP-DFTB method requested"
    end select

    status = int(dftb_status, c_int64_t)
    if (status /= 0_c_int64_t) then
      call write_c_string(message, status_message, status_message_len)
      return
    endif

    reference_energy = result%total_energy
    if (.not. compute_gradient) then
      call write_c_string(message, status_message, status_message_len)
      return
    endif

    if (.not. allocated(dftb_gradient)) then
      status = -9004_c_int64_t
      call write_c_string("OpenQP-DFTB gradient was requested but not returned", &
          status_message, status_message_len)
      return
    endif
    do iatom = 1, natom
      gradient(:, iatom) = real(dftb_gradient(:, iatom), c_double)
    enddo
    call write_c_string(message, status_message, status_message_len)
  end subroutine run_openqp_dftb

  subroutine set_spin_flip_reference(options)
    type(openqp_dftb_options_t), intent(inout) :: options

    options%scc = .true.
    options%spin_polarized = .true.
    options%spin_complete = .false.
    options%mrsf = .false.
    options%multiplicity = 3_dftb_i8
  end subroutine set_spin_flip_reference

  subroutine check_response_root(root, n_roots, status)
    integer(c_int64_t), intent(in) :: root
    integer(dftb_i8), intent(in) :: n_roots
    integer(dftb_i8), intent(inout) :: status

    if (root < 1_c_int64_t .or. root > int(n_roots, c_int64_t)) status = -9005_dftb_i8
  end subroutine check_response_root

  subroutine collect_unique_atomic_numbers(atoms, unique_atomic_numbers)
    type(openqp_dftb_atom_t), intent(in) :: atoms(:)
    integer(dftb_i8), allocatable, intent(out) :: unique_atomic_numbers(:)
    integer(dftb_i8), allocatable :: scratch(:)
    integer :: iatom
    integer :: n_unique

    allocate(scratch(size(atoms)))
    n_unique = 0
    do iatom = 1, size(atoms)
      if (.not. any(scratch(1:n_unique) == atoms(iatom)%atomic_number)) then
        n_unique = n_unique + 1
        scratch(n_unique) = atoms(iatom)%atomic_number
      endif
    enddo
    allocate(unique_atomic_numbers(n_unique))
    unique_atomic_numbers = scratch(1:n_unique)
  end subroutine collect_unique_atomic_numbers

  function c_string(chars, n_chars) result(text)
    character(kind=c_char), intent(in) :: chars(*)
    integer(c_int32_t), value :: n_chars
    character(len=:), allocatable :: text
    integer :: i
    integer :: n

    n = max(0, int(n_chars))
    allocate(character(len=n) :: text)
    do i = 1, n
      text(i:i) = chars(i)
    enddo
  end function c_string

  function lowercase(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i
    integer :: code

    lower = text
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar("A") .and. code <= iachar("Z")) lower(i:i) = achar(code + 32)
    enddo
  end function lowercase
#endif

end module openqp_dftb_bridge_mod
