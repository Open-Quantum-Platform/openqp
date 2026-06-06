"""SCF iteration table must print pathological values as numbers, not '*****'.

The table previously used bare f17.10/f14.8 edit descriptors; any value too
large for the column (e.g. the garbage energies produced by the dangling
ecp_zn buffer) printed as asterisks, hiding exactly the information needed to
debug.  The columns now route through fmt_real17/fmt_real14, which keep the
usual fixed-point form (and column width) while the value fits and fall back
to same-width scientific notation otherwise.
"""

import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_iteration_table_routes_reals_through_width_preserving_formatters():
    source = (ROOT / "source" / "scf.F90").read_text()

    assert "function fmt_real17" in source
    assert "function fmt_real14" in source
    assert "fmt_real17(energy%etot)" in source
    assert "fmt_real14(diis_error)" in source


def test_formatters_fall_back_to_same_width_scientific_notation():
    source = (ROOT / "source" / "scf.F90").read_text()

    # fixed-point first choice, scientific fallback, identical column widths
    assert "(f17.10)" in source
    assert "(es17.8e3)" in source
    assert "(f14.8)" in source
    assert "(es14.5e3)" in source


def test_formatter_descriptors_print_numbers_not_asterisks(tmp_path):
    """Compile a tiny formatter check so overflow fallback is behavioral."""
    gfortran = shutil.which("gfortran") or shutil.which("gfortran-15")
    if not gfortran:
        pytest.skip("gfortran is not available")

    source = tmp_path / "check_scf_format.f90"
    exe = tmp_path / "check_scf_format"
    source.write_text(
        textwrap.dedent(
            """
            program check_scf_format
              implicit none
              integer, parameter :: dp = kind(1.0d0)
              character(len=17) :: e_fixed, e_large
              character(len=14) :: g_fixed, g_large

              e_fixed = fmt_real17(-74.9631468000_dp)
              e_large = fmt_real17(-315184587.3_dp)
              g_fixed = fmt_real14(0.12345678_dp)
              g_large = fmt_real14(123456.0_dp)

              if (len(e_fixed) /= 17 .or. len(e_large) /= 17) error stop 1
              if (len(g_fixed) /= 14 .or. len(g_large) /= 14) error stop 2
              if (index(e_fixed, '*') /= 0 .or. index(e_large, '*') /= 0) error stop 3
              if (index(g_fixed, '*') /= 0 .or. index(g_large, '*') /= 0) error stop 4
              if (index(e_fixed, 'E') /= 0 .or. index(g_fixed, 'E') /= 0) error stop 5
              if (index(e_large, 'E') == 0 .or. index(g_large, 'E') == 0) error stop 6
            contains
              function fmt_real17(val) result(str)
                real(kind=dp), intent(in) :: val
                character(len=17) :: str
                if (val == val .and. abs(val) < 1.0e5_dp) then
                  write(str, '(f17.10)') val
                else
                  write(str, '(es17.8e3)') val
                end if
              end function fmt_real17

              function fmt_real14(val) result(str)
                real(kind=dp), intent(in) :: val
                character(len=14) :: str
                if (val == val .and. abs(val) < 1.0e4_dp) then
                  write(str, '(f14.8)') val
                else
                  write(str, '(es14.5e3)') val
                end if
              end function fmt_real14
            end program check_scf_format
            """
        )
    )

    subprocess.run([gfortran, str(source), "-o", str(exe)], check=True)
    subprocess.run([str(exe)], check=True)
