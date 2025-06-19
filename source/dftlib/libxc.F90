!> @brief  MODULE libxc
!> @brief  The head of libxc driver
!> @author Igor S. Gerasimov
!> @date   July, 2019 - Initial release -
!> @date   July, 2021 Making internal subroutines private
!> @todo   add LC-, CAM- free coefficient functionals
!> @todo   add meta-GGA functionals with laplacian of electron density
module libxc
  use xc_f03_lib_m
  use xc_f03_funcs_m
  use functionals, only: functional_t
  use precision, only: fp
  implicit none
  private
  character(len=*), parameter :: ABORTING = "Aborting from LibXC interface..."
  public :: libxc_input, libxc_destroy
contains

  !> @brief  setting up of using libxc functionals
  !>         (Analog of INPGDFT)
  !> @detail
  !          For adding a new functional:
  !           1) search code of interested functional that is a part of LibXC (like PBEX has code 101 and PBEC has code 130)
  !           2) Determine coefficients of each needed coefficients (as for example PBE-0.1: 0.1HF + 0.9PBEX + 1.0PBEC)
  !           3) Create a name of your functional (like PBEH). It can has any symbol excepting newline and space.
  !           4) Do not forget to set flag needtau=.true. if your functional need it.
  !
  !           the example of code:
  !             case("PBEH", "PBE-0.1") !two variant of names
  !                 HFEX = 0.10_fp !0.1HF exchange
  !                 call functional%add_functional(XC_GGA_X_PBE,0.9_fp) !PBE_X*0.9
  !                 call functional%add_functional(XC_GGA_C_PBE,1.0_fp) !PBE_C*1.0
  !> @author Igor S. Gerasimov
  !> @date   July, 2019 - Initial release -
  !> @date   July, 2021 Using messages module
  !>                    Adding optional arguments
  !> @date   Dec,  2022 Pass functional instead of using global
  !> @params functional_name (in) functional's name
  !> @param  infos           (inout)  info datatype
  !> @param  functional      (inout)  constructing functional
  subroutine libxc_input(functional_name, dft_params, tddft_params, functional)
    use messages, only: show_message, WITH_ABORT
    use types, only: dft_parameters, tddft_parameters

    character(len=*), intent(in) :: functional_name
    type(dft_parameters), intent(inout) :: dft_params
    type(tddft_parameters), intent(inout) :: tddft_params
    type(functional_t), intent(inout) :: functional

    ! Functional names
    character(len=:), allocatable :: funcname

    ! LibXC strings
    character(len=1024) :: LibXC_DOI, LibXC_reference, LibXC_version

    real(kind=fp) :: HFEX, MP2
    real(kind=fp) :: chf, c2opp, c2same
    HFEX = 0.0_fp
    chf = 0.0_fp
    c2opp = 0.0_fp
    c2same = 0.0_fp

    call xc_f03_reference(LibXC_reference)
    call xc_f03_reference_doi(LibXC_DOI)
    call xc_f03_version_string(LibXC_version)

    call show_message(" ") !empty line
    call show_message(" The LibXC "//trim(LibXC_version)//" version is used.")
    call show_message(" "//trim(LibXC_reference))
    call show_message(" The libXC interfaces are described in the following article:")
    call show_message(" Igor S. Gerasimov, Federico Zahariev, Sarom S. Leang, Anton Tesliuk, Mark S. Gordon, Michael G. Medvedev,")
    call show_message(" Introducing LibXC into GAMESS (US),")
    call show_message(" Mendeleev Commun., 2021, 31, 302–305")
    call show_message(" ") !empty line
    call show_message(" The information about selected functionals:")

    funcname = functional_name
    select case (funcname)
    case ("HFEX")
      HFEX = 1.00_fp
    case ("SLATER")
      call functional%add_functional(XC_LDA_X, 1.00_fp) !SLATER_X
    case ("TETER")
      call functional%add_functional(XC_LDA_XC_TETER93, 1.00_fp) !TETER93_XC
    case ("KSDT")
      call functional%add_functional(XC_LDA_XC_KSDT, 1.00_fp) !KSDT_XC
    case ("CORRKSDT")
      call functional%add_functional(XC_LDA_XC_CORRKSDT, 1.00_fp)
    case ("GDSMFB")
      call functional%add_functional(XC_LDA_XC_GDSMFB, 1.00_fp) !GDSMFB_XC
    case ("ZLP-20")
      call functional%add_functional(XC_LDA_XC_ZLP, 1.00_fp) !ZLP_XC (LDA)
    case ("PBE", "PBEPBE")
      call functional%add_functional(XC_GGA_X_PBE, 1.00_fp) !PBE_X
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp) !PBE_C
    case ("REVPBE")
      call functional%add_functional(XC_GGA_X_PBE_R, 1.00_fp) !PBE_R_X
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp) !PBE_C
    case ("RPBE")
      call functional%add_functional(XC_GGA_X_RPBE, 1.00_fp) !RPBE_X
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp) !PBE_C
    case ("PW91")
      call functional%add_functional(XC_GGA_X_PW91, 1.00_fp)
      call functional%add_functional(XC_GGA_C_PW91, 1.00_fp)
    case ("AM05")
      call functional%add_functional(XC_GGA_X_AM05, 1.00_fp)
      call functional%add_functional(XC_GGA_C_AM05, 1.00_fp)
    case ("PBESOL")
      call functional%add_functional(XC_GGA_X_PBE_SOL, 1.00_fp) !PBEsol_X
      call functional%add_functional(XC_GGA_C_PBE_SOL, 1.00_fp) !PBEsol_C
    case ("WC")
      call functional%add_functional(XC_GGA_X_WC, 1.00_fp)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp)
    case ("CHACHIYO")
      call functional%add_functional(XC_GGA_X_CHACHIYO, 1.00_fp)
      call functional%add_functional(XC_GGA_C_CHACHIYO, 1.00_fp)
    case ("BLYP")
      call functional%add_functional(XC_GGA_X_B88, 1.00_fp)
      call functional%add_functional(XC_GGA_C_LYP, 1.00_fp)
    case ("OLYP")
      call functional%add_functional(XC_GGA_X_OPTX, 1.00_fp)
      call functional%add_functional(XC_GGA_C_LYP, 1.00_fp)
    case ("XLYP")
      call functional%add_functional(XC_GGA_XC_XLYP, 1.00_fp)
    case ("KT1")
      call functional%add_functional(XC_GGA_XC_KT1, 1.00_fp)
    case ("KT2")
      call functional%add_functional(XC_GGA_XC_KT2, 1.00_fp)
    case ("KT3")
      call functional%add_functional(XC_GGA_XC_KT3, 1.00_fp)
    case ("TH1")
      call functional%add_functional(XC_GGA_XC_TH1, 1.00_fp)
    case ("TH2")
      call functional%add_functional(XC_GGA_XC_TH2, 1.00_fp)
    case ("TH3")
      call functional%add_functional(XC_GGA_XC_TH3, 1.00_fp)
    case ("TH4")
      call functional%add_functional(XC_GGA_XC_TH4, 1.00_fp)
    case ("HCTH93")
      call functional%add_functional(XC_GGA_XC_HCTH_93, 1.00_fp)
    case ("HCTH120")
      call functional%add_functional(XC_GGA_XC_HCTH_120, 1.00_fp)
    case ("HCTH147")
      call functional%add_functional(XC_GGA_XC_HCTH_147, 1.00_fp)
    case ("HCTH407")
      call functional%add_functional(XC_GGA_XC_HCTH_407, 1.00_fp)
    case ("HCTH407+", "HCTH407P")
      call functional%add_functional(XC_GGA_XC_HCTH_407P, 1.00_fp)
    case ("HCTH-A")
      call functional%add_functional(XC_GGA_X_HCTH_A, 1.00_fp)
      call functional%add_functional(XC_GGA_C_HCTH_A, 1.00_fp)
    case ("HCTH-P14")
      call functional%add_functional(XC_GGA_XC_HCTH_P14, 1.00_fp)
    case ("HCTH-P76")
      call functional%add_functional(XC_GGA_XC_HCTH_P76, 1.00_fp)
    case ("PBE1W")
      call functional%add_functional(XC_GGA_XC_PBE1W, 1.00_fp)
    case ("MPWLYP1W")
      call functional%add_functional(XC_GGA_XC_MPWLYP1W, 1.00_fp)
    case ("PBELYP1W")
      call functional%add_functional(XC_GGA_XC_PBELYP1W, 1.00_fp)
    case ("GAM")
      call functional%add_functional(XC_GGA_X_GAM, 1.00_fp)
      call functional%add_functional(XC_GGA_C_GAM, 1.00_fp)
    case ("MOHLYP")
      call functional%add_functional(XC_GGA_XC_MOHLYP, 1.00_fp)
    case ("MOHLYP-2", "MOHLYP2")
      call functional%add_functional(XC_GGA_XC_MOHLYP2, 1.00_fp)
    case ("SOGGA11")
      call functional%add_functional(XC_GGA_X_SOGGA11, 1.00_fp)
      call functional%add_functional(XC_GGA_C_SOGGA11, 1.00_fp)
    case ("N12")
      call functional%add_functional(XC_GGA_C_N12, 1.00_fp)
      call functional%add_functional(XC_GGA_X_N12, 1.00_fp)
    case ("HLE16")
      call functional%add_functional(XC_GGA_XC_HLE16, 1.00_fp)
    case ("EDF1")
      call functional%add_functional(XC_GGA_XC_EDF1, 1.00_fp)
    case ("LB07")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_LB07, 1.00_fp,&
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, omega=dft_params%cam_mu)
    case ("NCAP")
      call functional%add_functional(XC_GGA_XC_NCAP, 1.00_fp)
    case ("HTBS")
      call functional%add_functional(XC_GGA_X_HTBS, 1.00_fp)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp)
    case ("PKZB")
      call functional%add_functional(XC_MGGA_X_PKZB, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_PKZB, 1.00_fp)
    case ("TPSS")
      call functional%add_functional(XC_MGGA_X_TPSS, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_TPSS, 1.00_fp)
    case ("REVTPSS")
      call functional%add_functional(XC_MGGA_X_REVTPSS, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_REVTPSS, 1.00_fp)
    case ("MODTPSS")
      call functional%add_functional(XC_MGGA_X_MODTPSS, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_TPSS, 1.00_fp)
    case ("TPSSLOC", "TPSS-LOC")
      call functional%add_functional(XC_MGGA_X_TPSS, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_TPSSLOC, 1.00_fp)
    case ("RTPSS")
      call functional%add_functional(XC_MGGA_X_RTPSS, 1.00_fp)
      call functional%add_functional(XC_GGA_C_REGTPSS, 1.00_fp)
    case ("REGTPSS")
      call functional%add_functional(XC_MGGA_X_REGTPSS, 1.00_fp)
      call functional%add_functional(XC_GGA_C_REGTPSS, 1.00_fp)
    case ("MVS")
      call functional%add_functional(XC_MGGA_X_MVS, 1.00_fp)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp)
    case ("MVSB")
      call functional%add_functional(XC_MGGA_X_MVSB, 1.00_fp)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp)
    case ("MVSB*", "MVSBS")
      call functional%add_functional(XC_MGGA_X_MVSBS, 1.00_fp)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp)
    case ("MS0")
      call functional%add_functional(XC_MGGA_X_MS0, 1.00_fp)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp)
    case ("MS1")
      call functional%add_functional(XC_MGGA_X_MS1, 1.00_fp)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp)
    case ("MS2")
      call functional%add_functional(XC_MGGA_X_MS2, 1.00_fp)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp)
    case ("MS2B")
      call functional%add_functional(XC_MGGA_X_MS2B, 1.00_fp)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp)
    case ("MS2B*", "MS2BS")
      call functional%add_functional(XC_MGGA_X_MS2BS, 1.00_fp)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp)
    case ("SCAN")
      call functional%add_functional(XC_MGGA_X_SCAN, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_SCAN, 1.00_fp)
    case ("RSCAN", "REGSCAN")
      call functional%add_functional(XC_MGGA_X_RSCAN, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_RSCAN, 1.00_fp)
    case ("RPPSCAN", "R++SCAN")
      call functional%add_functional(XC_MGGA_X_RPPSCAN, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_RPPSCAN, 1.00_fp)
    case ("R2SCAN")
      call functional%add_functional(XC_MGGA_X_R2SCAN, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_R2SCAN, 1.00_fp)
    case ("R2SCAN01")
      call functional%add_functional(XC_MGGA_X_R2SCAN01, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_R2SCAN01, 1.00_fp)
    case ("R4SCAN")
      call functional%add_functional(XC_MGGA_X_R4SCAN, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_R2SCAN, 1.00_fp)
    case ("REVSCAN")
      call functional%add_functional(XC_MGGA_X_REVSCAN, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_REVSCAN, 1.00_fp)
    case ("TASK")
      call functional%add_functional(XC_MGGA_X_TASK, 1.00_fp)
      call functional%add_functional(XC_LDA_C_PW, 1.00_fp)
    case ("TM")
      call functional%add_functional(XC_MGGA_X_TM, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_TM, 1.00_fp)
    case ("REVTM")
      call functional%add_functional(XC_MGGA_X_REVTM, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_REVTM, 1.00_fp)
    case ("REGTM")
      call functional%add_functional(XC_MGGA_X_REGTM, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_TM, 1.00_fp)
    case ("RREGTM")
      call functional%add_functional(XC_MGGA_X_REGTM, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_RREGTM, 1.00_fp)
    case ("MGGAC")
      call functional%add_functional(XC_MGGA_X_MGGAC, 1.00_fp)
      call functional%add_functional(XC_GGA_C_MGGAC, 1.00_fp)
    case ("RMGGAC")
      call functional%add_functional(XC_MGGA_X_MGGAC, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_RMGGAC, 1.00_fp)
    case ("TAUHCTH", "THCTH")
      call functional%add_functional(XC_MGGA_X_TAU_HCTH, 1.00_fp)
      call functional%add_functional(XC_GGA_C_TAU_HCTH, 1.00_fp)
    case ("VSXC")
      call functional%add_functional(XC_MGGA_X_GVT4, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_VSXC, 1.00_fp)
    case ("TPSSLYP1W")
      call functional%add_functional(XC_MGGA_XC_TPSSLYP1W, 1.00_fp)
    case ("M06-L", "M06L")
      call functional%add_functional(XC_MGGA_X_M06_L, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_M06_L, 1.00_fp)
    case ("REVM06-L", "REVM06L")
      call functional%add_functional(XC_MGGA_X_REVM06_L, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_REVM06_L, 1.00_fp)
    case ("M11-L", "M11L")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_MGGA_X_M11_L, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
      call functional%add_functional(XC_MGGA_C_M11_L, 1.00_fp)
    case ("MN12-L", "MN12L")
      call functional%add_functional(XC_MGGA_X_MN12_L, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_MN12_L, 1.00_fp)
    case ("MN15-L", "MN15L")
      call functional%add_functional(XC_MGGA_X_MN15_L, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_MN15_L, 1.00_fp)
    case ("HLE17")
      call functional%add_functional(XC_MGGA_XC_HLE17, 1.00_fp)
    case ("HLTA")
      call functional%add_functional(XC_MGGA_X_HLTA, 1.00_fp)
      call functional%add_functional(XC_MGGA_C_HLTAPW, 1.00_fp)
    case ("LDA0")
      call functional%add_functional(XC_HYB_LDA_XC_LDA0, 1.00_fp, hfex=HFEX) !LDA0_hXC
    case ("CAM-LDA0") ! The energy of He in 3-21G basis set looks OK (-2.8529482409)
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_LDA_XC_CAM_LDA0, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu) !CAM-LDA0_hXC
    case ("APF")
      call functional%add_functional(XC_HYB_GGA_XC_APF, 1.00_fp, hfex=HFEX)
    case ("B1WC")
      call functional%add_functional(XC_HYB_GGA_XC_B1WC, 1.00_fp, hfex=HFEX)
    case ("PBE0", "PBE-25", "PBE1PBE")
      call functional%add_functional(XC_HYB_GGA_XC_PBEH, 1.00_fp, hfex=HFEX)
    case ("PBE-33")
      call functional%add_functional(XC_HYB_GGA_XC_PBE0_13, 1.00_fp, hfex=HFEX)
    case ("PBE-38", "PBE-3/8")
      call functional%add_functional(XC_HYB_GGA_XC_PBE38, 1.00_fp, hfex=HFEX)
    case ("PBE-50")
      call functional%add_functional(XC_HYB_GGA_XC_PBE50, 1.00_fp, hfex=HFEX)
    case ("PBE-2X")
      call functional%add_functional(XC_HYB_GGA_XC_PBE_2X, 1.00_fp, hfex=HFEX)
    case ("PBE-3X")
      HFEX = 0.84_fp
      call functional%add_functional(XC_GGA_X_PBE, 0.16_fp) !PBE_X
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp) !PBE_C
    case ("BHH")
      call functional%add_functional(XC_HYB_GGA_XC_BHANDH, 1.00_fp, hfex=HFEX)
    case ("B1PW91")
      call functional%add_functional(XC_HYB_GGA_XC_B1PW91, 1.00_fp, hfex=HFEX)
    case ("B3PW91")
      call functional%add_functional(XC_HYB_GGA_XC_B3PW91, 1.00_fp, hfex=HFEX)
    case ("BLYP35")
      call functional%add_functional(XC_HYB_GGA_XC_BLYP35, 1.00_fp, hfex=HFEX)
    case ("B1LYP")
      call functional%add_functional(XC_HYB_GGA_XC_B1LYP, 1.00_fp, hfex=HFEX)
    case ("BHHLYP")
      call functional%add_functional(XC_HYB_GGA_XC_BHANDHLYP, 1.00_fp, hfex=HFEX)
    case ("STG1X")
      call functional%add_functional(XC_HYB_GGA_XC_BHANDHLYP, 1.00_fp, hfex=HFEX) 
      hfex=0.85
      tddft_params%spc_coco = 0.5_fp
      tddft_params%spc_ovov = 0.5_fp
      tddft_params%spc_coov = 0.5_fp
      write(*,fmt='(a)') "STG1X = B(0.15)HF(0.85)-LYP functional"
      write(*,fmt='(3a)') "[3] Y. Horbatenko, S. Lee, M. Filatov, and C. H. Choi, ", &
            "J. Phys. Chem. A, 123, 7991-8000 (2019); ", &
            "DOI: 10.1021/acs.jpca.9b07556"
    case ("B3LYP")
      call show_message("B3LYP functional has different meanings,")
      call show_message("so that it can not be run using LibXC interface.")
      call show_message("For running B3LYP, choose one of them:")
      call show_message(" - B3LYPV1R with  VWN RPA LDA correlation part (default for Gaussian)")
      call show_message(" - B3LYPV3  with  VWN_3   LDA correlation part")
      call show_message(" - B3LYPV5  with  VWN_5   LDA correlation part (default for GAMESS-US)")
      call show_message(ABORTING, WITH_ABORT)
    case ("B3LYPV1R")
      call functional%add_functional(XC_HYB_GGA_XC_B3LYP, 1.00_fp, hfex=HFEX)
    case ("B3LYPV3")
      call functional%add_functional(XC_HYB_GGA_XC_B3LYP3, 1.00_fp, hfex=HFEX)
    case ("B3LYPV5", "B3LYP5")
      call functional%add_functional(XC_HYB_GGA_XC_B3LYP5, 1.00_fp, hfex=HFEX)
    case ("O3LYP")
      call functional%add_functional(XC_HYB_GGA_XC_O3LYP, 1.00_fp, hfex=HFEX)
    case ("X3LYP")
      call functional%add_functional(XC_HYB_GGA_XC_X3LYP, 1.00_fp, hfex=HFEX)
    case ("REVB3LYP")
      call functional%add_functional(XC_HYB_GGA_XC_REVB3LYP, 1.00_fp, hfex=HFEX)
    case ("B3LYP*", "B3LYPS")
      call functional%add_functional(XC_HYB_GGA_XC_B3LYPs, 1.00_fp, hfex=HFEX)
    case ("B50LYP", "B5050LYP")
      call functional%add_functional(XC_HYB_GGA_XC_B5050LYP, 1.00_fp, hfex=HFEX)
    case ("KMLYP")
      call functional%add_functional(XC_HYB_GGA_XC_KMLYP, 1.00_fp, hfex=HFEX)
    case ("CASE21")
      call functional%add_functional(XC_HYB_GGA_XC_CASE21, 1.00_fp, hfex=HFEX)
    case ("B97-1P")
      call functional%add_functional(XC_HYB_GGA_XC_B97_1p, 1.00_fp, hfex=HFEX)
    case ("B97")
      call functional%add_functional(XC_HYB_GGA_XC_B97, 1.00_fp, hfex=HFEX)
    case ("B97-1")
      call functional%add_functional(XC_HYB_GGA_XC_B97_1, 1.00_fp, hfex=HFEX)
    case ("B97-2")
      call functional%add_functional(XC_HYB_GGA_XC_B97_2, 1.00_fp, hfex=HFEX)
    case ("B97-3")
      call functional%add_functional(XC_HYB_GGA_XC_B97_3, 1.00_fp, hfex=HFEX)
    case ("B97-K")
      call functional%add_functional(XC_HYB_GGA_XC_B97_K, 1.00_fp, hfex=HFEX)
    case ("MPWLYP1M")
      call functional%add_functional(XC_HYB_GGA_XC_MPWLYP1M, 1.00_fp, hfex=HFEX)
    case ("MPW1K")
      call functional%add_functional(XC_HYB_GGA_XC_MPW1K, 1.00_fp, hfex=HFEX)
    case ("MPW1PBE")
      call functional%add_functional(XC_HYB_GGA_XC_MPW1PBE, 1.00_fp, hfex=HFEX)
    case ("MPW1PW", "MPW1PW91")
      call functional%add_functional(XC_HYB_GGA_XC_MPW1PW, 1.00_fp, hfex=HFEX)
    case ("MPW1LYP")
      call functional%add_functional(XC_HYB_GGA_XC_MPW1LYP, 1.00_fp, hfex=HFEX)
    case ("MPW3PW", "MPW3PW91")
      call functional%add_functional(XC_HYB_GGA_XC_MPW3PW, 1.00_fp, hfex=HFEX)
    case ("MPW3LYP")
      call functional%add_functional(XC_HYB_GGA_XC_MPW3LYP, 1.00_fp, hfex=HFEX)
    case ("SOGGA11-X", "SOGGA11X")
      call functional%add_functional(XC_HYB_GGA_X_SOGGA11_X, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_GGA_C_SOGGA11_X, 1.00_fp)
    case ("QTP17")
      call functional%add_functional(XC_HYB_GGA_XC_QTP17, 1.00_fp, hfex=HFEX)
    case ("EDF2")
      call functional%add_functional(XC_HYB_GGA_XC_EDF2, 1.00_fp, hfex=HFEX)
    case ("HPBEINT")
      call functional%add_functional(XC_HYB_GGA_XC_HPBEINT, 1.00_fp, hfex=HFEX)
    case ("CAP0")
      call functional%add_functional(XC_HYB_GGA_XC_CAP0, 1.00_fp, hfex=HFEX)
    case ("WC04")
      call functional%add_functional(XC_HYB_GGA_XC_WC04, 1.00_fp, hfex=HFEX)
    case ("WP04")
      call functional%add_functional(XC_HYB_GGA_XC_WP04, 1.00_fp, hfex=HFEX)
    case ("WB97")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_WB97, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("WB97X")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_WB97X, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("WB97X-D")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_WB97X_D, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("HSE03")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_HSE03, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("HSE06")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_HSE06, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("HSE12")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_HSE12, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("HSE12-S")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_HSE12S, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("HSESOL")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_HSE_SOL, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("MCAM-B3LYP", "MCAMB3LYP")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_MCAM_B3LYP, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("CAM-B3LYP", "CAMB3LYP")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_CAM_B3LYP, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("CAMH-B3LYP", "CAMHB3LYP")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_CAMH_B3LYP, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("DTCAM-TUNE", "CDTCAMTUNE")  ! Use to tune HF exchange in DFT and TDDFT
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_TUNED_CAM_B3LYP, 1.00_fp, &
            external_parameters=(/ 0.81_fp, &
                                   dft_params%cam_alpha+dft_params%cam_beta, &
                                  -dft_params%cam_beta, &
                                   dft_params%cam_mu /), &
            alpha=dft_params%cam_alpha, &
            beta=dft_params%cam_beta, &
            omega=dft_params%cam_mu)
      write(*,fmt='(3a)') "[2] W. Park, A. Lashkaripour, K. Komarov, S. Lee, M. Huix-Rotllant, ", &
            "and C. H. Choi, J. Chem. Theory Comput., 20(13), 5679-5694 (2024); ", &
            "DOI: 10.1021/acs.jctc.4c00640"
    case ("DTCAM-VAEE", "DTCAMVAEE")  ! see doi.org/10.1021/acs.jctc.4c00640
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_TUNED_CAM_B3LYP, 1.00_fp, &
            external_parameters=(/ 0.81_fp, 0.30_fp,  0.20_fp, 0.33_fp  /), &
            alpha=dft_params%cam_alpha, & ! = 0.50
            beta=dft_params%cam_beta, &   ! =-0.20
            omega=dft_params%cam_mu)      ! = 0.33
      tddft_params%cam_alpha = 0.5_fp
      tddft_params%cam_beta = -0.1_fp
      tddft_params%cam_mu = dft_params%cam_mu
      tddft_params%spc_coco = 0.5_fp
      tddft_params%spc_ovov = 0.5_fp
      tddft_params%spc_coov = 0.5_fp
      write(*,fmt='(3a)') "[2] W. Park, A. Lashkaripour, K. Komarov, S. Lee, M. Huix-Rotllant, ", &
            "and C. H. Choi, J. Chem. Theory Comput., 20(13), 5679-5694 (2024); ", &
            "DOI: 10.1021/acs.jctc.4c00640"
    case ("DTCAM-XIV", "DTCAMXIV")  ! see doi.org/10.1021/acs.jctc.4c00640
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_TUNED_CAM_B3LYP, 1.00_fp, &
            external_parameters=(/ 0.81_fp, 0.30_fp,  0.29_fp, 0.33_fp  /), &
            alpha=dft_params%cam_alpha, & ! = 0.59
            beta=dft_params%cam_beta, &   ! =-0.29
            omega=dft_params%cam_mu)      ! = 0.33
      tddft_params%cam_alpha = 0.50_fp
      tddft_params%cam_alpha = 0.10_fp
      tddft_params%cam_beta = 0.90_fp
      tddft_params%cam_mu = dft_params%cam_mu
      tddft_params%spc_coco = 0.50_fp
      tddft_params%spc_ovov = 0.50_fp
      tddft_params%spc_coov = 0.50_fp
      write(*,fmt='(3a)') "[2] W. Park, A. Lashkaripour, K. Komarov, S. Lee, M. Huix-Rotllant, ", &
            "and C. H. Choi, J. Chem. Theory Comput., 20(13), 5679-5694 (2024); ", &
            "DOI: 10.1021/acs.jctc.4c00640"
    case ("DTCAM-XI", "DTCAMXI")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_TUNED_CAM_B3LYP, 1.00_fp, &
            external_parameters=(/ 0.81_fp, 1.02_fp, -0.52_fp, 0.33_fp   /), &
            alpha=dft_params%cam_alpha, & ! = 0.50
            beta=dft_params%cam_beta, &   ! = 0.52
            omega=dft_params%cam_mu)      ! = 0.33
      tddft_params%cam_alpha = 0.395_fp
      tddft_params%cam_beta = 0.425_fp
      tddft_params%cam_mu = dft_params%cam_mu
      tddft_params%spc_coco = 0.50_fp
      tddft_params%spc_ovov = 0.50_fp
      tddft_params%spc_coov = 0.50_fp
      write(*,fmt='(3a)') "[2] W. Park, A. Lashkaripour, K. Komarov, S. Lee, M. Huix-Rotllant, ", &
            "and C. H. Choi, J. Chem. Theory Comput., ??, ?? (2024); ", &
            "DOI: 10.1021/acs.jctc.4c00640"
    case ("DTCAM-AEE", "DTCAMAEE")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_TUNED_CAM_B3LYP, 1.00_fp, &
            external_parameters=(/ 0.81_fp, 0.48_fp, -0.29_fp, 0.33_fp  /), &
            alpha=dft_params%cam_alpha, & ! = 0.19
            beta=dft_params%cam_beta, &   ! = 0.29
            omega=dft_params%cam_mu)      ! = 0.33
      tddft_params%cam_alpha = 0.15_fp
      tddft_params%cam_beta = 0.95_fp
      tddft_params%cam_mu = dft_params%cam_mu
      tddft_params%spc_coco = 0.50_fp
      tddft_params%spc_ovov = 0.50_fp
      tddft_params%spc_coov = 0.50_fp
      write(*,fmt='(3a)') "[2] K. Komarov, W. Park, S. Lee, M. Huix-Rotllant, ", &
            "and C. H. Choi, J. Chem. Theory Comput., 19, 7671-7684 (2023); ", &
            "DOI: 10.1021/acs.jctc.3c00884"
    case ("DTCAM-VEE", "DTCAMVEE")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_TUNED_CAM_B3LYP, 1.00_fp, &
            external_parameters=(/ 0.81_fp, 0.48_fp, -0.29_fp, 0.33_fp  /), &
            alpha=dft_params%cam_alpha, & ! = 0.19
            beta=dft_params%cam_beta, &   ! = 0.29
            omega=dft_params%cam_mu)      ! = 0.33
      tddft_params%cam_alpha = 0.48_fp
      tddft_params%cam_beta =  0.00_fp
      tddft_params%cam_mu = dft_params%cam_mu
      tddft_params%spc_coco = 0.50_fp
      tddft_params%spc_ovov = 0.50_fp
      tddft_params%spc_coov = 0.50_fp
      write(*,fmt='(3a)') "[2] K. Komarov, W. Park, S. Lee, M. Huix-Rotllant, ", &
            "and C. H. Choi, J. Chem. Theory Comput., 19, 7671-7684 (2023); ", &
            "DOI: 10.1021/acs.jctc.3c00884"
    case ("DTCAM-STG", "DTCAMSTG")  
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_TUNED_CAM_B3LYP, 1.00_fp, &
            external_parameters=(/ 0.81_fp, 0.24_fp,  0.24_fp, 0.33_fp  /), &
            alpha=dft_params%cam_alpha, & ! = 0.48
            beta=dft_params%cam_beta, &   ! =-0.24
            omega=dft_params%cam_mu)      ! = 0.33
      tddft_params%cam_alpha = 0.71_fp
      tddft_params%cam_beta = -0.1_fp
      tddft_params%cam_mu = dft_params%cam_mu
      tddft_params%spc_coco = 0.5_fp
      tddft_params%spc_ovov = 0.5_fp
      tddft_params%spc_coov = 0.5_fp
      write(*,fmt='(3a)') "[2] W. Park, A. Lashkaripour, K. Komarov, S. Lee, M. Huix-Rotllant, ", &
            "and C. H. Choi, J. Chem. Theory Comput., 20(13), 5679-5694 (2024); ", &
            "DOI: 10.1021/acs.jctc.4c00640"
    case ("RCAM-B3LYP", "RCAMB3LYP")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_RCAM_B3LYP, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("TUNCAM-B3LYP", "TUNEDCAM-B3LYP")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_TUNED_CAM_B3LYP, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("CAM-QTP00", "CAMQTP00", "CAM-QTP(00)", "CAM-QTP-00")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_CAM_QTP_00, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("CAM-QTP01", "CAMQTP01", "CAM-QTP(01)", "CAM-QTP-01")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_CAM_QTP_01, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("CAM-QTP02", "CAMQTP02", "CAM-QTP(02)", "CAM-QTP-02")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_CAM_QTP_02, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("CAM-PBEH") ! The energy of He in 3-21G basis set looks OK (-2.8610170571)
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_CAM_PBEH, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("LRC-WPBEH")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_LRC_WPBEH, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("LRC-WPBE")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_LRC_WPBE, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("LC-WPBE", "LCWPBE")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_LC_WPBE, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("WHPBE0")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_WHPBE0, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("N12-SX", "N12SX")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_X_N12_SX, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
      call functional%add_functional(XC_GGA_C_N12_SX, 1.00_fp)
    case ("LC-QTP")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_LC_QTP, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("LC-BLYP")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_GGA_XC_LC_BLYP, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
    case ("TPSSH")
      call functional%add_functional(XC_HYB_MGGA_XC_TPSSH, 1.00_fp, hfex=HFEX)
    case ("TPSS0")
      call functional%add_functional(XC_HYB_MGGA_XC_TPSS0, 1.00_fp, hfex=HFEX)
    case ("REVTPSSH")
      call functional%add_functional(XC_HYB_MGGA_XC_REVTPSSH, 1.00_fp, hfex=HFEX)
    case ("MS2H")
      call functional%add_functional(XC_HYB_MGGA_X_MS2H, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp)
    case ("MVSH")
      call functional%add_functional(XC_HYB_MGGA_X_MVSH, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp)
    case ("SCAN0")
      call functional%add_functional(XC_HYB_MGGA_X_SCAN0, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_MGGA_C_SCAN, 1.00_fp)
    case ("REVSCAN0")
      call functional%add_functional(XC_HYB_MGGA_X_REVSCAN0, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_MGGA_C_REVSCAN, 1.00_fp)
    case ("HYBTAUHCTH", "HYBTHCTH", "HYB-THCTH", "HYB-TAUHCTH", "THCTHHYB", "TAUHCTHHYB")
      call functional%add_functional(XC_HYB_MGGA_X_TAU_HCTH, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_GGA_C_HYB_TAU_HCTH, 1.00_fp)
    case ("BMK")
      call functional%add_functional(XC_HYB_MGGA_X_BMK, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_GGA_C_BMK, 1.00_fp)
    case ("DLDF")
      call functional%add_functional(XC_HYB_MGGA_X_DLDF, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_MGGA_C_DLDF, 1.00_fp)
    case ("M05")
      call functional%add_functional(XC_HYB_MGGA_X_M05, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_MGGA_C_M05, 1.00_fp)
    case ("M05-2X", "M052X")
      call functional%add_functional(XC_HYB_MGGA_X_M05_2X, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_MGGA_C_M05_2X, 1.00_fp)
    case ("M06")
      call functional%add_functional(XC_HYB_MGGA_X_M06, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_MGGA_C_M06, 1.00_fp)
    case ("REVM06")
      call functional%add_functional(XC_HYB_MGGA_X_REVM06, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_MGGA_C_REVM06, 1.00_fp)
    case ("M06-2X", "M062X")
      call functional%add_functional(XC_HYB_MGGA_X_M06_2X, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_MGGA_C_M06_2X, 1.00_fp)
    case ("M06-HF", "M06HF")
      call functional%add_functional(XC_HYB_MGGA_X_M06_HF, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_MGGA_C_M06_HF, 1.00_fp)
    case ("M08-SO")
      call functional%add_functional(XC_HYB_MGGA_X_M08_SO, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_MGGA_C_M08_SO, 1.00_fp)
    case ("M08-HX")
      call functional%add_functional(XC_HYB_MGGA_X_M08_HX, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_MGGA_C_M08_HX, 1.00_fp)
    case ("MN15")
      call functional%add_functional(XC_HYB_MGGA_X_MN15, 1.00_fp, hfex=HFEX)
      call functional%add_functional(XC_MGGA_C_MN15, 1.00_fp)
    case ("PW86BC95", "PW86B95")
      call functional%add_functional(XC_HYB_MGGA_XC_PW86B95, 1.00_fp, hfex=HFEX)
    case ("BB1K")
      call functional%add_functional(XC_HYB_MGGA_XC_BB1K, 1.00_fp, hfex=HFEX)
    case ("MPW1BC95", "MPW1B95")
      call functional%add_functional(XC_HYB_MGGA_XC_MPW1B95, 1.00_fp, hfex=HFEX)
    case ("MPW1B1K")
      call functional%add_functional(XC_HYB_MGGA_XC_MPWB1K, 1.00_fp, hfex=HFEX)
    case ("PW6BC95", "PW6B95")
      call functional%add_functional(XC_HYB_MGGA_XC_PW6B95, 1.00_fp, hfex=HFEX)
    case ("PWB6K")
      call functional%add_functional(XC_HYB_MGGA_XC_PWB6K, 1.00_fp, hfex=HFEX)
    case ("MPW1KCIS")
      call functional%add_functional(XC_HYB_MGGA_XC_MPW1KCIS, 1.00_fp, hfex=HFEX)
    case ("MPWKCIS1K")
      call functional%add_functional(XC_HYB_MGGA_XC_MPWKCIS1K, 1.00_fp, hfex=HFEX)
    case ("B0KCIS")
      call functional%add_functional(XC_HYB_MGGA_XC_B0KCIS, 1.00_fp, hfex=HFEX)
    case ("PBE1KCIS")
      call functional%add_functional(XC_HYB_MGGA_XC_PBE1KCIS, 1.00_fp, hfex=HFEX)
    case ("TPSS1KCIS")
      call functional%add_functional(XC_HYB_MGGA_XC_TPSS1KCIS, 1.00_fp, hfex=HFEX)
    case ("X1BC95", "X1B95")
      call functional%add_functional(XC_HYB_MGGA_XC_X1B95, 1.00_fp, hfex=HFEX)
    case ("XB1K")
      call functional%add_functional(XC_HYB_MGGA_XC_XB1K, 1.00_fp, hfex=HFEX)
    case ("B88BC95", "BBC95", "B88B95", "BB95")
      call functional%add_functional(XC_HYB_MGGA_XC_B88B95, 1.00_fp, hfex=HFEX)
    case ("B86BC95", "B86B95")
      call functional%add_functional(XC_HYB_MGGA_XC_B86B95, 1.00_fp, hfex=HFEX)
    case ("M06-SX", "M06SX")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_MGGA_X_M06_SX, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
      call functional%add_functional(XC_MGGA_C_M06_SX, 1.00_fp)
    case ("M11")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_MGGA_X_M11, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
      call functional%add_functional(XC_MGGA_C_M11, 1.00_fp)
    case ("REVM11")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_MGGA_X_REVM11, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
      call functional%add_functional(XC_MGGA_C_REVM11, 1.00_fp)
    case ("MN12-SX", "MN12SX")
      dft_params%cam_flag = .true.
      call functional%add_functional(XC_HYB_MGGA_X_MN12_SX, 1.00_fp, &
              alpha=dft_params%cam_alpha, &
              beta=dft_params%cam_beta, &
              omega=dft_params%cam_mu)
      call functional%add_functional(XC_MGGA_C_MN12_SX, 1.00_fp)
    case ("PBE0-DH")
      dft_params%dh_flag = .true.
      HFEX = 0.5_fp
      MP2 = 0.128_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_GGA_X_PBE, 1.00_fp-HFEX)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp-MP2)
    case ("TPSS0-DH")
      dft_params%dh_flag = .true.
      HFEX = 0.5_fp
      MP2 = 0.128_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_MGGA_X_TPSS, 1.00_fp-HFEX)
      call functional%add_functional(XC_MGGA_C_TPSS, 1.00_fp-MP2)
    case ("SCAN0-DH")
      dft_params%dh_flag = .true.
      HFEX = 0.5_fp
      MP2 = 0.128_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_MGGA_X_SCAN, 1.00_fp-HFEX)
      call functional%add_functional(XC_MGGA_C_SCAN, 1.00_fp-MP2)
    case ("PBE-QIDH")
      dft_params%dh_flag = .true.
      HFEX = 0.693391_fp
      MP2 = 0.128_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_GGA_X_PBE, 1.00_fp-HFEX)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp-MP2)
    case ("TPSS-QIDH")
      dft_params%dh_flag = .true.
      HFEX = 0.693391_fp
      MP2 = 0.128_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_MGGA_X_TPSS, 1.00_fp-HFEX)
      call functional%add_functional(XC_MGGA_C_TPSS, 1.00_fp-MP2)
    case ("SCAN-QIDH")
      dft_params%dh_flag = .true.
      HFEX = 0.693391_fp
      MP2 = 0.128_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_MGGA_X_SCAN, 1.00_fp-HFEX)
      call functional%add_functional(XC_MGGA_C_SCAN, 1.00_fp-MP2)
    case ("PBE0-2")
      dft_params%dh_flag = .true.
      HFEX = 0.793701_fp
      MP2 = 0.5_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_GGA_X_PBE, 1.00_fp-HFEX)
      call functional%add_functional(XC_GGA_C_PBE, 1.00_fp-MP2)
    case ("TPSS0-2")
      dft_params%dh_flag = .true.
      HFEX = 0.793701_fp
      MP2 = 0.5_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_MGGA_X_TPSS, 1.00_fp-HFEX)
      call functional%add_functional(XC_MGGA_C_TPSS, 1.00_fp-MP2)
    case ("SCAN0-2")
      dft_params%dh_flag = .true.
      HFEX = 0.793701_fp
      MP2 = 0.5_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_MGGA_X_SCAN, 1.00_fp-HFEX)
      call functional%add_functional(XC_MGGA_C_SCAN, 1.00_fp-MP2)
    case ("B2-PLYP", "B2PLYP")
      dft_params%dh_flag = .true.
      HFEX = 0.53_fp
      MP2 = 0.27_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_GGA_X_B88, 1.00_fp-HFEX) !Becke88
      call functional%add_functional(XC_GGA_C_LYP, 1.00_fp-MP2) !Lee-Yang-Parr
    case ("B2GP-PLYP", "B2GPPLYP")
      dft_params%dh_flag = .true.
      HFEX = 0.65_fp
      MP2 = 0.36_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_GGA_X_B88, 1.00_fp-HFEX) !Becke88
      call functional%add_functional(XC_GGA_C_LYP, 1.00_fp-MP2) !Lee-Yang-Parr
    case ("B2T-PLYP", "B2TPLYP")
      dft_params%dh_flag = .true.
      HFEX = 0.60_fp
      MP2 = 0.31_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_GGA_X_B88, 1.00_fp-HFEX) !Becke88
      call functional%add_functional(XC_GGA_C_LYP, 1.00_fp-MP2) !Lee-Yang-Parr
    case ("B2K-PLYP", "B2KPLYP")
      dft_params%dh_flag = .true.
      HFEX = 0.42_fp
      MP2 = 0.72_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_GGA_X_B88, 1.00_fp-HFEX) !Becke88
      call functional%add_functional(XC_GGA_C_LYP, 1.00_fp-MP2) !Lee-Yang-Parr
    case ("MPW2-PLYP", "MPW2PLYP")
      dft_params%dh_flag = .true.
      HFEX = 0.53_fp
      MP2 = 0.27_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_GGA_X_MPW91, 1.00_fp-HFEX) !mPW91
      call functional%add_functional(XC_GGA_C_LYP, 1.00_fp-MP2) !Lee-Yang-Parr
    case ("MPWK-PLYP", "MPWKPLYP")
      dft_params%dh_flag = .true.
      HFEX = 0.42_fp
      MP2 = 0.72_fp
      chf = HFEX
      c2same = MP2
      c2opp = MP2
      call functional%add_functional(XC_GGA_X_MPW91, 1.00_fp-HFEX) !mPW91
      call functional%add_functional(XC_GGA_C_LYP, 1.00_fp-MP2) !Lee-Yang-Parr
    case default
      !It will be good, if someone write procedure for parsing of functional names, like PBE-1-PBE or PBE-0 or S+PBE-3-LYP+VWN-RPA
      call show_message("Unrecognized functional name: "//funcname)
      call show_message("Please, check the documentation about this functional or")
      call show_message(" implement it to LibXC [https://gitlab.com/libxc/libxc],")
      call show_message("   and, then, to OQP-LibXC interface")
      call show_message(ABORTING, WITH_ABORT)
    end select
    ! No one of functionals was not selected
    if (.not. functional%can_calculate() .and. funcname /= "HFEX") then
      call show_message("No one functional was not selected!")
      call show_message("Please, check the input!")
      call show_message(ABORTING, WITH_ABORT)
    end if
    dft_params%HFScale = HFEX
    dft_params%MP2SS_Scale = c2opp
    dft_params%MP2OS_Scale = c2same
    ! Print information about non-local part of functionals
    call show_message(" ")
    if (.not. dft_params%cam_flag) then
      call show_message("(A,ES16.8E2)", "The global hybrid part:                ", dft_params%HFScale)
    else
      call show_message("CAM-corrected functional is called")
      call show_message("Yanai's notation of CAM parameters is used")
      call show_message("(A,ES16.8E2)", "The global hybrid part:                ", dft_params%cam_alpha)
      call show_message("(A,ES16.8E2)", "Additional long-range hybrid part:     ", dft_params%cam_beta)
      call show_message("(A,ES16.8E2)", "Range-saparated factor:                ", dft_params%cam_mu)
    end if
    if (dft_params%dh_flag) then
      call show_message(" ")
      call show_message("Double-hybrid functional is called")
      call show_message("(A,ES16.8E2)", "Same-spin MP2 correlation factor:      ", dft_params%MP2SS_Scale)
      call show_message("(A,ES16.8E2)", "Opposite-spin MP2 correlation factor:  ", dft_params%MP2OS_Scale)
    end if
    call show_message(" ")
  end subroutine libxc_input
  !
  !> @brief  Destroy internal variables of functional
  !> @author Igor S. Gerasimov
  !> @date   Dec, 2020 - Initial release -
  !> @date   Dec, 2022 Pass local functional instead of global
  subroutine libxc_destroy(functional)
    type(functional_t), intent(inout) :: functional
    call functional%destroy
  end subroutine libxc_destroy
end module libxc
