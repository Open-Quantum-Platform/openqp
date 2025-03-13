#include "tagarray.fh"
module oqp_tagarray_driver
  use tagarray

  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_char, c_ptr, c_null_ptr

  implicit none
  private
  character(len=*), parameter, private :: module_name = "oqp_tagarray_driver"
  public :: tagarray_get_cptr
  character(len=*), parameter, public :: OQP_prefix = "OQP::"
  character(len=*), parameter, public :: OQP_DM_A = OQP_prefix // "DM_A"
  character(len=*), parameter, public :: OQP_DM_B = OQP_prefix // "DM_B"
  character(len=*), parameter, public :: OQP_FOCK_A = OQP_prefix // "FOCK_A"
  character(len=*), parameter, public :: OQP_FOCK_B = OQP_prefix // "FOCK_B"
  character(len=*), parameter, public :: OQP_E_MO_A = OQP_prefix // "E_MO_A"
  character(len=*), parameter, public :: OQP_E_MO_B = OQP_prefix // "E_MO_B"
  character(len=*), parameter, public :: OQP_VEC_MO_A = OQP_prefix // "VEC_MO_A"
  character(len=*), parameter, public :: OQP_VEC_MO_B = OQP_prefix // "VEC_MO_B"
  character(len=*), parameter, public :: OQP_Hcore = OQP_prefix // "Hcore"
  character(len=*), parameter, public :: OQP_SM = OQP_prefix // "SM"
  character(len=*), parameter, public :: OQP_TM = OQP_prefix // "TM"
  character(len=*), parameter, public :: OQP_WAO = OQP_prefix // "WAO"
  character(len=*), parameter, public :: OQP_td_abxc = OQP_prefix // "td_abxc"
  character(len=*), parameter, public :: OQP_td_bvec_mo = OQP_prefix // "td_bvec_mo"
  character(len=*), parameter, public :: OQP_td_mrsf_density = OQP_prefix // "td_mrsf_density"
  character(len=*), parameter, public :: OQP_td_p = OQP_prefix // "td_p"
  character(len=*), parameter, public :: OQP_td_t = OQP_prefix // "td_t"
  character(len=*), parameter, public :: OQP_td_xpy = OQP_prefix // "td_xpy"
  character(len=*), parameter, public :: OQP_td_xmy = OQP_prefix // "td_xmy"
  character(len=*), parameter, public :: OQP_td_energies = OQP_prefix // "td_energies"
  character(len=*), parameter, public :: OQP_log_filename = OQP_prefix // "log_filename"
  character(len=*), parameter, public :: OQP_basis_filename = OQP_prefix // "basis_filename"
  character(len=*), parameter, public :: OQP_hbasis_filename = OQP_prefix // "hbasis_filename"
! Used to compute properties between two geometries
  character(len=*), parameter, public :: OQP_xyz_old = OQP_prefix // "xyz_old"
  character(len=*), parameter, public :: OQP_overlap_ao = OQP_prefix // "overlap_ao_non_orthogonal"
  character(len=*), parameter, public :: OQP_overlap_mo = OQP_prefix // "overlap_mo_non_orthogonal"
  character(len=*), parameter, public :: OQP_E_MO_A_old = OQP_prefix // "E_MO_A_old"
  character(len=*), parameter, public :: OQP_E_MO_B_old = OQP_prefix // "E_MO_B_old"
  character(len=*), parameter, public :: OQP_VEC_MO_A_old = OQP_prefix // "VEC_MO_A_old"
  character(len=*), parameter, public :: OQP_VEC_MO_B_old = OQP_prefix // "VEC_MO_B_old"
  character(len=*), parameter, public :: OQP_td_bvec_mo_old = OQP_prefix // "td_bvec_mo_old"
  character(len=*), parameter, public :: OQP_td_energies_old = OQP_prefix // "td_energies_old"
  character(len=*), parameter, public :: OQP_nac = OQP_prefix // "nac"
  character(len=*), parameter, public :: OQP_td_states_phase = OQP_prefix // "td_states_phase"
  character(len=*), parameter, public :: OQP_td_states_overlap = OQP_prefix // "td_states_overlap"

  character(len=*), parameter, public :: OQP_DM_A_comment = "Alpha-spin triangle Density matrix"
  character(len=*), parameter, public :: OQP_DM_B_comment = "Beta-spin triangle Density matrix"
  character(len=*), parameter, public :: OQP_FOCK_A_comment = "Alpha-spin triangle Fock matrix"
  character(len=*), parameter, public :: OQP_FOCK_B_comment = "Beta-spin triangle Fock matrix"
  character(len=*), parameter, public :: OQP_E_MO_A_comment = "Energies of alpha molecular orbitals"
  character(len=*), parameter, public :: OQP_E_MO_B_comment = "Energies of beta molecular orbitals"
  character(len=*), parameter, public :: OQP_VEC_MO_A_comment = "Coefficients of alpha molecular orbitals"
  character(len=*), parameter, public :: OQP_VEC_MO_B_comment = "Coefficients of beta molecular orbitals"
  character(len=*), parameter, public :: OQP_Hcore_comment = "triangle core Hamiltonian matrix"
  character(len=*), parameter, public :: OQP_SM_comment = "triangle Overlap matrix"
  character(len=*), parameter, public :: OQP_TM_comment = "triangle Kinetic-Energy matrix"
  character(len=*), parameter, public :: OQP_WAO_comment = "??? WAO ???"
  character(len=*), parameter, public :: OQP_td_abxc_comment = "??? td_abxc ???"
  character(len=*), parameter, public :: OQP_td_bvec_mo_comment = "??? td_bvec_mo ???"
  character(len=*), parameter, public :: OQP_td_mrsf_density_comment = "??? td_mrsf_density ???"
  character(len=*), parameter, public :: OQP_td_p_comment = "??? td_p ???"
  character(len=*), parameter, public :: OQP_td_t_comment = "??? td_t ???"
  character(len=*), parameter, public :: OQP_td_xpy_comment = OQP_prefix // "(X+Y) vector for target state in TD-DFT calculations"
  character(len=*), parameter, public :: OQP_td_xmy_comment = OQP_prefix // "(X-Y) vector for target state in TD-DFT calculations"
  character(len=*), parameter, public :: OQP_td_energies_comment = OQP_prefix // "Responce energies"
  character(len=*), parameter, public :: OQP_log_filename_comment = OQP_prefix // "log filename"
  character(len=*), parameter, public :: OQP_basis_filename_comment = OQP_prefix // "basis filename"
  character(len=*), parameter, public :: OQP_hbasis_filename_comment = OQP_prefix // "Huckel basis_filename for Huckel Guess"
  character(len=*), parameter, public :: OQP_nac_comment = OQP_prefix // "nonadiabatic coupling nstates x nstates"
  character(len=*), parameter, public :: OQP_overlap_mo_comment = OQP_prefix // "overlap between MOs of geo1 and geo2"
  character(len=*), parameter, public :: OQP_overlap_ao_comment = OQP_prefix // "overlap between geo1 and geo2"
  character(len=*), parameter, public :: OQP_td_states_phase_comment = OQP_prefix // "Bvecs phase sign with respect to Bvec_old"
  character(len=*), parameter, public :: OQP_td_states_overlap_comment = OQP_prefix // "Bvecs phase sign with respect to Bvec_old"
  character(len=*), parameter, public :: OQP_xyz_oldcomment = OQP_prefix // "saved geo from previous step"
  character(len=*), parameter, public :: all_tags(32) = (/ character(len=80) :: &
    OQP_DM_A, OQP_DM_B, OQP_FOCK_A, OQP_FOCK_B, OQP_E_MO_A, OQP_E_MO_B, &
    OQP_VEC_MO_A, OQP_VEC_MO_B, OQP_Hcore, OQP_SM, OQP_TM, OQP_WAO, &
    OQP_td_abxc, OQP_td_bvec_mo, OQP_td_mrsf_density, OQP_td_p, OQP_td_t, &
    OQP_log_filename, OQP_basis_filename, OQP_hbasis_filename, &
    OQP_xyz_old, OQP_overlap_mo, OQP_overlap_ao, OQP_E_MO_A_old, OQP_E_MO_B_old, &
    OQP_VEC_MO_A_old, OQP_VEC_MO_B_old, OQP_td_bvec_mo_old, OQP_td_energies_old, &
    OQP_nac, OQP_td_states_phase, OQP_td_states_overlap /)
  interface tagarray_get_data
    module procedure tagarray_get_data_int64_val, tagarray_get_data_int64_1d, tagarray_get_data_int64_2d, tagarray_get_data_int64_3d
    module procedure tagarray_get_data_real64_val, tagarray_get_data_real64_1d, tagarray_get_data_real64_2d, tagarray_get_data_real64_3d
    module procedure tagarray_get_data_char8_val, tagarray_get_data_char8_1d
  end interface
  interface data_has_tags
    module procedure data_has_tags_location, data_has_tags_ms
  end interface
  public :: data_has_tags, check_status
  public :: tagarray_get_data
  public :: TA_TYPE_INT64, TA_TYPE_REAL64, TA_TYPE_CHAR8
  public :: ta_ok
contains

    function tagarray_get_cptr(container, tag, ptr, type_id, ndims, dims, data_size) result(res)
    type(container_t), intent(inout) :: container
    character(len=*), intent(in) :: tag
    type(c_ptr), intent(out) :: ptr
    integer(c_int64_t) :: res
    integer(c_int32_t), optional, intent(out) :: type_id
    integer(c_int32_t), optional, intent(out) :: ndims
    integer(c_int64_t), optional, intent(out) :: dims(:)
    integer(c_int64_t), optional, intent(out) :: data_size

    type(recordinfo_t) :: record_info

    ptr = c_null_ptr
    record_info = container%get_record_info(tag)
    res = container%get_status()

    if (res == TA_OK) then
      ptr = record_info%data
      res = product(record_info%dimensions(1:record_info%n_dimensions))
      if (present(type_id)) type_id = record_info%type_id
      if (present(ndims  )) ndims   = record_info%n_dimensions
      if (present(dims   )) dims    = record_info%dimensions
      if (present(data_size   )) data_size    = record_info%data_length
    end if

  end function tagarray_get_cptr


  subroutine data_has_tags_location(container, tags, location, abort, status)
    use messages, only: show_message, WITHOUT_ABORT
    type(container_t), intent(inout) :: container
    character(len=*), intent(in) :: tags(:)
    character(len=*), intent(in) :: location
    logical, optional, intent(in) :: abort
    integer(c_int32_t), optional, intent(out) :: status
    integer(c_int32_t) :: tag_id, status_
    logical :: abort_
    abort_ = WITHOUT_ABORT
    if (present(abort)) abort_ = abort
    status_ = container%has_records(tags, tag_id)
    if (status_ /= TA_OK) call show_message( &
        location // ": " // get_status_message(status_, trim(tags(tag_id))), &
        abort_)
    if (present(status)) status = status_
  end subroutine data_has_tags_location
  subroutine data_has_tags_ms(container, tags, modulename, subroutinename, abort, status)
    use messages, only: show_message, WITHOUT_ABORT
    type(container_t), intent(inout) :: container
    character(len=*), intent(in) :: tags(:)
    character(len=*), intent(in) :: modulename, subroutinename
    logical, optional, intent(in) :: abort
    integer(c_int32_t), optional, intent(out) :: status
    integer(c_int32_t) :: tag_id, status_
    logical :: abort_
    abort_ = WITHOUT_ABORT
    if (present(abort)) abort_ = abort
    status_ = container%has_records(tags, tag_id)
    if (status_ /= TA_OK) call show_message( &
        modulename // "::" // subroutinename // ": " // get_status_message(status_, trim(tags(tag_id))), &
        abort_)
    if (present(status)) status = status_
  end subroutine data_has_tags_ms
  subroutine check_status(status, modulename, subroutinename, tag, abort)
    use messages, only: show_message, WITHOUT_ABORT
    integer(c_int32_t), intent(in) :: status
    character(len=*), intent(in) :: modulename, subroutinename, tag
    logical, optional, intent(in) :: abort
    logical :: abort_
    abort_ = WITHOUT_ABORT
    if (present(abort)) abort_ = abort
    if (status /= TA_OK) call show_message( &
        modulename // "::" // subroutinename // ": " // get_status_message(status, trim(tag)), &
        abort_)
  end subroutine check_status
  subroutine tagarray_get_data_int64_val(container, tag, ptr, status)
    type(container_t), intent(inout) :: container
    character(len=*), intent(in) :: tag
    integer(8), pointer :: ptr
    integer(c_int32_t), optional, intent(out) :: status
    integer(c_int32_t) :: status_
    TA_GET_CONTAINER_VALUE(container, tag, TA_TYPE_INT64, ptr, status_)
    if (present(status)) status = status_
  end subroutine tagarray_get_data_int64_val
  subroutine tagarray_get_data_int64_1d(container, tag, ptr, status)
    type(container_t), intent(inout) :: container
    character(len=*), intent(in) :: tag
    integer(8), pointer :: ptr(:)
    integer(c_int32_t), optional, intent(out) :: status
    integer(c_int32_t) :: status_
    TA_GET_CONTAINER_DATA(container, tag, TA_TYPE_INT64, ptr, status_)
    if (present(status)) status = status_
  end subroutine tagarray_get_data_int64_1d
  subroutine tagarray_get_data_int64_2d(container, tag, ptr, status)
    type(container_t), intent(inout) :: container
    character(len=*), intent(in) :: tag
    integer(8), pointer :: ptr(:,:)
    integer(c_int32_t), optional, intent(out) :: status
    integer(c_int32_t) :: status_
    TA_GET_CONTAINER_DATA(container, tag, TA_TYPE_INT64, ptr, status_)
    if (present(status)) status = status_
  end subroutine tagarray_get_data_int64_2d
  subroutine tagarray_get_data_int64_3d(container, tag, ptr, status)
    type(container_t), intent(inout) :: container
    character(len=*), intent(in) :: tag
    integer(8), pointer :: ptr(:,:,:)
    integer(c_int32_t), optional, intent(out) :: status
    integer(c_int32_t) :: status_
    TA_GET_CONTAINER_DATA(container, tag, TA_TYPE_INT64, ptr, status_)
    if (present(status)) status = status_
  end subroutine tagarray_get_data_int64_3d
  subroutine tagarray_get_data_real64_val(container, tag, ptr, status)
    type(container_t), intent(inout) :: container
    character(len=*), intent(in) :: tag
    real(8), pointer :: ptr
    integer(c_int32_t), optional, intent(out) :: status
    integer(c_int32_t) :: status_
    TA_GET_CONTAINER_VALUE(container, tag, TA_TYPE_REAL64, ptr, status_)
    if (present(status)) status = status_
  end subroutine tagarray_get_data_real64_val
  subroutine tagarray_get_data_real64_1d(container, tag, ptr, status)
    type(container_t), intent(inout) :: container
    character(len=*), intent(in) :: tag
    real(8), pointer :: ptr(:)
    integer(c_int32_t), optional, intent(out) :: status
    integer(c_int32_t) :: status_
    TA_GET_CONTAINER_DATA(container, tag, TA_TYPE_REAL64, ptr, status_)
    if (present(status)) status = status_
  end subroutine tagarray_get_data_real64_1d
  subroutine tagarray_get_data_real64_2d(container, tag, ptr, status)
    type(container_t), intent(inout) :: container
    character(len=*), intent(in) :: tag
    real(8), pointer :: ptr(:,:)
    integer(c_int32_t), optional, intent(out) :: status
    integer(c_int32_t) :: status_
    TA_GET_CONTAINER_DATA(container, tag, TA_TYPE_REAL64, ptr, status_)
    if (present(status)) status = status_
  end subroutine tagarray_get_data_real64_2d
  subroutine tagarray_get_data_real64_3d(container, tag, ptr, status)
    type(container_t), intent(inout) :: container
    character(len=*), intent(in) :: tag
    real(8), pointer :: ptr(:,:,:)
    integer(c_int32_t), optional, intent(out) :: status
    integer(c_int32_t) :: status_
    TA_GET_CONTAINER_DATA(container, tag, TA_TYPE_REAL64, ptr, status_)
    if (present(status)) status = status_
  end subroutine tagarray_get_data_real64_3d

  subroutine tagarray_get_data_char8_val(container, tag, ptr, status)
    type(container_t), intent(inout) :: container
    character(len=*), intent(in) :: tag
    character(len=*, kind=c_char), pointer :: ptr
    integer(c_int32_t), optional, intent(out) :: status
    integer(c_int32_t) :: status_
    TA_GET_CONTAINER_VALUE(container, tag, TA_TYPE_CHAR8, ptr, status_)
    if (present(status)) status = status_
  end subroutine tagarray_get_data_char8_val

  subroutine tagarray_get_data_char8_1d(container, tag, ptr, status)
    type(container_t), intent(inout) :: container
    character(len=*), intent(in) :: tag
    character(len=*, kind=c_char), pointer :: ptr(:)
    integer(c_int32_t), optional, intent(out) :: status
    integer(c_int32_t) :: status_
    TA_GET_CONTAINER_DATA(container, tag, TA_TYPE_CHAR8, ptr, status_)
    if (present(status)) status = status_
  end subroutine tagarray_get_data_char8_1d
end module oqp_tagarray_driver
