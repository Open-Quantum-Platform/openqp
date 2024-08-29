#!/usr/bin/python

"""
>    This script generates 
>      1) libxc_choose_functional subroutine for libxc.src from all available LibXC functionals
>      2) documentation of libxc_choose_functional subroutine for docs-input.txt
"""

import re
import pylibxc

text = {
"LDA_K":"LDA kinetic",
"GGA_K":"GGA kinetic",
"MGGA_K":"meta-GGA kinetic",
"LDA_X":"LDA exchange",
"GGA_X":"GGA exchange",
"MGGA_X":"meta-GGA exchange",
"LDA_C":"LDA correlation",
"GGA_C":"GGA correlation",
"MGGA_C":"meta-GGA correlation",
"LDA_XC":"LDA exchange-correlation",
"GGA_XC":"GGA exchange-correlation",
"MGGA_XC":"meta-GGA exchange-correlation",
"HLDA_X":"hybrid LDA exchange",
"HGGA_X":"hybrid GGA exchange",
"HMGGA_X":"hybrid meta-GGA exchange",
"HLDA_XC":"hybrid LDA exchange-correlation",
"HGGA_XC":"hybrid GGA exchange-correlation",
"HMGGA_XC":"hybrid meta-GGA exchange-correlation"
}

sequence = ["LDA_K", "LDA_X", "LDA_C", "LDA_XC", "HLDA_X", "HLDA_XC",
            "GGA_K", "GGA_X", "GGA_C", "GGA_XC", "HGGA_X", "HGGA_XC",
            "MGGA_K","MGGA_X","MGGA_C","MGGA_XC","HMGGA_X","HMGGA_XC"]

xc_OK     = []
out       = []
variables = []
values    = []
code      = []
deinit    = []

xc_func_list = pylibxc.util.xc_available_functional_names()

for xc_func in xc_func_list:
  func = pylibxc.LibXCFunctional(xc_func, "unpolarized")
  knd = func.get_kind()
  fgs = func.get_flags()
  if knd > pylibxc.flags.XC_EXCHANGE_CORRELATION:
    print("{} was removed: it is not ex/cor/ex-cor functional".format(xc_func))
    continue
  if not fgs & pylibxc.flags.XC_FLAGS_3D:
    print("{} was removed: it is not 3D functional".format(xc_func))
    continue
  if fgs & pylibxc.flags.XC_FLAGS_NEEDS_LAPLACIAN:
    print("{} was removed: it needs Laplacian".format(xc_func))
    continue
  if fgs & pylibxc.flags.XC_FLAGS_VV10:
    print("{} was removed: it has VV10 correlation".format(xc_func))
    continue
  if fgs & pylibxc.flags.XC_FLAGS_HYB_LC:
    print("{} was removed: mixing of several LC/CAM functionals is not possible".format(xc_func))
    continue
  if fgs & pylibxc.flags.XC_FLAGS_HYB_CAM:
    print("{} was removed: mixing of several LC/CAM functionals is not possible".format(xc_func))
    continue
  if xc_func.upper() == "LDA_X_YUKAWA":
    print("{} was removed: it has Yukawa potential".format(xc_func))
    continue
  if fgs & pylibxc.flags.XC_FLAGS_HYB_LCY:
    print("{} was removed: it has Yukawa potential".format(xc_func))
    continue
  if fgs & pylibxc.flags.XC_FLAGS_HYB_CAMY:
    print("{} was removed: it has Yukawa potential".format(xc_func))
    continue
  if "VDW" in func.get_name().upper():
    print("{} was removed: GAMESS does not support vdW functionals".format(xc_func))
    continue
  xc_out = xc_func.upper().replace("HYB_", "H")
  if xc_out == "LDA_X":
    xc_out = "LDA_X_SLATER"
  elif xc_out == "LDA_C_XALPHA":
    xc_out = "LDA_X_XALPHA"
  xc_OK.append([xc_out, xc_func.upper(), func.get_name(), func.get_number()])

xc_OK.sort(key=lambda x: x[0])

func_types = list(set([ "_".join(x[0].split("_")[:2]) for x in xc_OK ]))
func_types.sort()

out.append("  subroutine libxc_choose_functionals(functional, read_unit)")
out.append("    use mod_nameio, only: input_group")
out.append("    type(functional_t), intent(inout) :: functional")
out.append("    integer, intent(in) :: read_unit")
out.append("    !")
out.append("    ! internal variables")

for i in func_types:
  fvarname = i + "_name"
  group = i + "_group"
  variables.append("    character(len={0:2}),  parameter   :: {1:30} = \"{2}\"".format(len(i), fvarname, i))
  variables.append("    class(input_group), allocatable :: {0}".format(group))
  code.append(     "    allocate({0})".format(group))
  code.append(     "    call {1}%read_section(read_unit, {0})".format(fvarname, group))
  deinit.append(   "    if(allocated({0})) then".format(group))
  deinit.append(   "      call {0}%destroy".format(group))
  deinit.append(   "      deallocate({0})".format(group))
  deinit.append(   "    end if")

for xc in xc_OK:
  func     = xc[0]
  funccode = "XC_" + xc[1]
  funcname = func + "_name"
  typ      = "_".join(func.split("_")[:2])
  typ_pair = typ + "_group"
  funcinp  = "_".join(func.split("_")[2:])
  variables.append("    character(len={0:2}),  parameter   :: {1:30} = {2:30} !{3:>4}: {4}".format(len(funcinp), funcname, '"' + funcinp.upper() + '"', xc[3], xc[2]))
  if "MGGA" not in typ:
    code.append( "    call add_free_functional(functional,{0},{1},{2})".format(typ_pair,funcname,funccode))
  else:
    code.append( "    call add_free_functional(functional,{0},{1},{2},.true.)".format(typ_pair,funcname,funccode))

out.extend(variables)
out.extend(values)
out.extend(code)
out.extend(deinit)
out.append("  end subroutine libxc_choose_functionals")

open("libxc_choose_functionals.src_part","w").write("\n".join(out))

out = []

for f in sequence:
  if(f[-1] == "K"): continue # skip kinetic functionals
  if(f not in func_types): continue
  out.append("")
  out.append("="*60)
  out.append("")
  out.append("${0:10} group        (relevant for DFTTYP=USELIBXC)".format(f))
  out.append("             (and also FUNCTIONAL is not set in $LIBXC)")
  out.append("")
  out.append("    This group sets {0}".format(text[f]))
  out.append("type of functionals with any coefficients as a double value.")
  out.append("Please, check the references of used functionals by yourself.")
  out.append("The available functionals in this group are next:")
  funcs = [ x for x in xc_OK if len(re.findall(r"^" + f + "_", x[0])) > 0 ]
  for func in funcs:
    out.append("        {0:20}= {1}".format("_".join(func[0].split("_")[2:]), func[2]))

out.append("")
open("docs-input.txt_part","w").write("\n".join(out))
