module soc_mrsf_mod

  implicit none

  character(len=*), parameter :: module_name = "soc_mrsf_mod"

  private

  public soc_mrsf

contains

  subroutine soc_mrsf_C(c_handle) bind(C, name="soc_mrsf")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call soc_mrsf(inf)
  end subroutine soc_mrsf_C

  subroutine soc_mrsf(infos)
    use io_constants, only: iw
    use types, only: information
    use oqp_tagarray_driver
    use precision, only: dp
    use physical_constants, only: UNITS_EV
    use printing, only: print_module_info
    use messages, only: show_message, with_abort

    implicit none

    character(len=*), parameter :: subroutine_name = "soc_mrsf"

    type(information), target, intent(inout) :: infos

    real(kind=dp), contiguous, pointer :: singlet_energies(:), triplet_energies(:)
    real(kind=dp), contiguous, pointer :: bvec_mo_s(:,:), bvec_mo_t(:,:)
    real(kind=dp) :: e_ref
    integer :: ns, nt, nstates, ist

    open(unit=iw, file=infos%log_filename, position="append")

    call print_module_info('SOC_MRSF', 'Spin-Orbit Coupling: MRSF Energies')

    e_ref = infos%mol_energy%energy

    call data_has_tags(infos%dat, &
        (/ character(len=80) :: OQP_td_singlet_energies, OQP_td_triplet_energies /), &
        module_name, subroutine_name, WITH_ABORT)

    call tagarray_get_data(infos%dat, OQP_td_singlet_energies, singlet_energies)
    call tagarray_get_data(infos%dat, OQP_td_triplet_energies, triplet_energies)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo_s, bvec_mo_s)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo_t, bvec_mo_t)

    ns = size(singlet_energies)
    nt = size(triplet_energies)
    nstates = max(ns, nt)

    write(iw,'(/,4x,a,f20.10,a)') 'Reference SCF energy: ', e_ref, ' Hartree'
    write(iw,'(/,4x,a5,2x,a20,2x,a20)') 'State', 'Singlet (au)', 'Triplet (au)'
    write(iw,'(4x,a5,2x,a20,2x,a20)') '-----', '--------------------', '--------------------'

    singlet_energies = singlet_energies + e_ref
    triplet_energies = triplet_energies + e_ref

    do ist = 1, nstates
      if (ist <= ns .and. ist <= nt) then
        write(iw,'(4x,i5,2x,f20.10,2x,f20.10)') ist, singlet_energies(ist), triplet_energies(ist)
      else if (ist <= ns) then
        write(iw,'(4x,i5,2x,f20.10,2x,a20)') ist, singlet_energies(ist), '---'
      else
        write(iw,'(4x,i5,2x,a20,2x,f20.10)') ist, '---', triplet_energies(ist)
      end if
    end do


    write(iw,'(a,2i5)') 'bvec_s shape:', size(bvec_mo_s,1), size(bvec_mo_s,2)
    write(iw,'(a,3f16.10)') 'bvec_s(1:3,1):', bvec_mo_s(1,1), bvec_mo_s(2,1), bvec_mo_s(3,1)
    write(iw,'(a,f16.10)') 'norm(bvec_s,1):', sqrt(dot_product(bvec_mo_s(:,1), bvec_mo_s(:,1)))

    write(iw,'(a,2i5)') 'bvec_s shape:', size(bvec_mo_t,1), size(bvec_mo_t,2)
    write(iw,'(a,3f16.10)') 'bvec_s(1:3,1):', bvec_mo_t(1,1), bvec_mo_t(2,1), bvec_mo_t(3,1)
    write(iw,'(a,f16.10)') 'norm(bvec_s,1):', sqrt(dot_product(bvec_mo_t(:,1), bvec_mo_t(:,1)))

    write(iw,'()')
    call flush(iw)
    close(iw)

  end subroutine soc_mrsf

end module soc_mrsf_mod
