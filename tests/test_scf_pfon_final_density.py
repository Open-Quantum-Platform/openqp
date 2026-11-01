"""A pFON convergence claim must use a density annealed to the final 1 K.

The temperature can drop from hundreds of kelvin to 1 K at the beginning of
an iteration.  Its energy and DIIS error still use the density built at the
previous temperature; accepting convergence on that iteration reports a hot
density as a 1 K solution.

Two decks exercise the guard:

* H2O ROHF/PBE/6-31G (examples/SCF/h2o_rohf_6-31g_pbe_pfon without
  save_mol).  It exposed the defect on Windows MKL ILP64 (0.00007740 Eh away
  from its reference at four OpenMP threads); on Linux MKL and macOS
  Accelerate its convergence does not coincide with the 1 K transition.
* Butadiene triplet ROHF/PBE/6-31G*.  On Linux MKL ILP64 the unfixed driver
  dropped from 1500 K to 1 K at the start of iteration 13 and exited there,
  reporting -155.6691111340 Eh from the 1500 K density.  The 1 K solution,
  reached after further cold iterations, is -155.6688516530 Eh (2.6e-4 Eh
  higher); this deck fails the temperature and the energy assertion on the
  unfixed code.
"""

import os
import re
import subprocess
import sys

import pytest


H2O = ("roks/pbe/6-31g mult=3 scf(maxit=200,pfon=true) dftgrid(pruned=none)\n",
       """\
3
H2O
O 0.0 0.0 -0.041061554
H -0.533194329 0.533194329 -0.614469223
H 0.533194329 -0.533194329 -0.614469223
""", None)

BUTADIENE = ("roks/pbe/6-31g* mult=3 scf(maxit=200,pfon=true)\n",
             """\
10
C4H6
C -1.901080641 0.114577198 0.0
C -0.574900598 -0.402237631 0.0
C 1.901080641 -0.114577198 0.0
C 0.574900598 0.402237631 0.0
H -2.755248674 -0.54511591 0.0
H -2.071562576 1.18308684 0.0
H 2.755248674 0.54511591 0.0
H -0.441348756 -1.477941036 0.0
H 2.071562576 -1.18308684 0.0
H 0.441348756 1.477941036 0.0
""", -155.6688516530)

# Well below the 2.6e-4 Eh error of the early exit and well above the
# 1e-8 Eh agreement of the 1 K solution across thread counts.
ENERGY_TOL = 1.0e-6


@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("case", [H2O, BUTADIENE], ids=["h2o", "butadiene"])
def test_pfon_energy_uses_final_temperature_density(tmp_path, threads, case):
    deck, geometry, reference = case
    (tmp_path / "mol.xyz").write_text(geometry)
    inp = tmp_path / "pfon.oqp"
    inp.write_text(deck + 'geom="mol.xyz"\n')
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    proc = subprocess.run(
        [sys.executable, "-c", "import sys; from oqp.pyoqp import main; sys.exit(main())",
         inp.name],
        cwd=tmp_path, env=env, text=True, capture_output=True, timeout=900)
    log_path = tmp_path / "pfon.log"
    log = log_path.read_text() if log_path.exists() else ""
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    assert "SCF convergence achieved" in log, log[-3000:]
    # The printed temperature is the one the iteration's new density is built
    # at; the Fock and energy of an iteration use the previous row's density.
    # Convergence may therefore be accepted only when the preceding row was
    # already at 1 K.
    temperatures = [float(x) for x in re.findall(r"Temp:\s+([0-9.]+)", log)]
    assert len(temperatures) >= 2, log[-3000:]
    assert temperatures[-2:] == [1.0, 1.0], temperatures[-8:]
    if reference is not None:
        energy = float(re.findall(r"TOTAL energy =\s+(-?\d+\.\d+)", log)[-1])
        assert abs(energy - reference) < ENERGY_TOL, (energy, reference)
