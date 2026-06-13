"""PyOQP

This module serves as the main entry point for the OQP (Open Quantum Package)
software. It handles command-line arguments, initializes the OQP environment,
and executes the appropriate computations or tests.
"""
import os
import sys
import time
import argparse
from signal import signal, SIGINT, SIG_DFL


def _apply_omp_threads_from_input(argv):
    """Honour an OpenMP thread-count request from the input file or CLI BEFORE
    the native OpenMP runtime initialises (it caches OMP_NUM_THREADS when liboqp
    loads, so this must run before `import oqp`).

    Sources, highest precedence first:
      * CLI:   --omp N   (or --omp=N)
      * input: a line  `omp_threads = N`  (typically in the [input] section)

    The value sets OMP_NUM_THREADS (threads per process / MPI rank).  If neither
    is given, the existing environment / built-in default is left untouched.
    """
    import re
    n = None
    for i, a in enumerate(argv):
        if a in ("--omp", "--omp-threads") and i + 1 < len(argv):
            n = argv[i + 1]
            break
        if a.startswith("--omp="):
            n = a.split("=", 1)[1]
            break
    if n is None:
        inp = next((a for a in argv[1:]
                    if not a.startswith("-") and os.path.isfile(a)), None)
        if inp:
            try:
                with open(inp, encoding="utf-8", errors="ignore") as fh:
                    m = re.search(r"(?mi)^[ \t]*omp_threads[ \t]*=[ \t]*(\d+)",
                                  fh.read())
                if m:
                    n = m.group(1)
            except OSError:
                pass
    if n is not None:
        try:
            ni = int(n)
        except (TypeError, ValueError):
            return
        if ni >= 1:
            os.environ["OMP_NUM_THREADS"] = str(ni)


def _set_threading_defaults():
    """Set conservative nested-threading defaults before native libraries load.

    OpenQP parallelizes the native integral and response kernels with OpenMP.
    BLAS must NOT spawn a second worker-thread layer inside those OMP regions,
    or it oversubscribes the cores (OMP_threads x BLAS_threads) and stalls --
    measured 1.53x faster on caffeine/6-31G(d,p) at 24 threads with sequential
    BLAS (31.1s -> 20.3s).  GNU libgomp worker stacks also need headroom for the
    integral kernels on macOS/arm64 (otherwise a startup segfault).  All values
    are setdefault, so an explicit user environment still wins.
    """
    defaults = {
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "OMP_STACKSIZE": "512M",
        "GOMP_STACKSIZE": "512M",
    }
    if sys.platform == "darwin":
        # macOS default: use the PERFORMANCE-core count, not all logical cores.
        # On Apple Silicon spilling onto the efficiency cores oversubscribes and
        # slows the integral build; on Intel Macs all cores are equal.  Derive it
        # per-host (hw.perflevel0.physicalcpu = P-cores on Apple Silicon, falling
        # back to hw.physicalcpu / os.cpu_count()) instead of a fixed cap, so we
        # never over- or under-subscribe a given machine.  setdefault, so an
        # explicit OMP_NUM_THREADS still wins.
        import subprocess
        ncore = None
        for key in ("hw.perflevel0.physicalcpu", "hw.physicalcpu"):
            try:
                ncore = int(subprocess.check_output(
                    ["sysctl", "-n", key], stderr=subprocess.DEVNULL).strip())
            except Exception:
                ncore = None
            if ncore and ncore > 0:
                break
        if not ncore or ncore < 1:
            ncore = os.cpu_count() or 1
        defaults["OMP_NUM_THREADS"] = str(ncore)

    for key, value in defaults.items():
        os.environ.setdefault(key, value)


# Establish conservative BLAS/OMP defaults first (setdefault, so they never
# override an explicit environment), then honour an explicit per-rank thread
# request from --omp / omp_threads, which hard-sets OMP_NUM_THREADS and wins.
_set_threading_defaults()
_apply_omp_threads_from_input(sys.argv)

import oqp
from oqp.utils.file_utils import dump_log
from oqp.utils.input_checker import check_input_values
from oqp.molecule import Molecule
from oqp.library.runfunc import (
   compute_energy, compute_grad, compute_nac, compute_soc, compute_geom, compute_md,
   compute_nacme, compute_properties, compute_data, compute_hess, compute_thermo,
   compute_namd
)
from oqp.utils.mpi_utils import MPIManager


class Runner:
    """
    OQP main class for running calculations and tests.

    one can read input parameters from text file using input_file, e.g,
        OQP(project=project_name, log=log_file, input_file=a_text_file)

    or

    one can pass input parameters to OQP internally using input_dict, e.g.,
        OQP(project=project_name, log=log_file, input_dict=USER_CONFIG)
        USER_CONFIG follows the same format as OQP_CONFIG_SCHEMA

    """

    def __init__(self, project=None, input_file=None, log=None,
                 input_dict=None, silent=0, usempi=True):
        """
        Initialize the OQP Runner.

        Args:
            project (str): Project name.
            input_file (str): Path to the input file.
            log (str): Path to the log file.
            input_dict (dict): Dictionary containing input parameters.
            silent (int): Flag to run in silent mode (0 or 1).
            usempi (bool): Flag to enable MPI functions
        """
        start_time = time.time()

        # Define the mapping of run types to their respective functions
        self.run_func = {
            'energy': compute_energy,
            'ekt': compute_energy,
            'grad': compute_grad,
            'nac': compute_nac,
            'nacme': compute_nacme,
            'soc': compute_soc,
            'optimize': compute_geom,
            'meci': compute_geom,
            'mecp': compute_geom,
            'mep': compute_geom,
            'md': compute_md,
            'namd': compute_namd,
            'ts': compute_geom,
            'tci': compute_geom,
            'irc': compute_geom,
            'neb': compute_geom,
            'hess': compute_hess,
            'thermo': compute_thermo,
            'prop': compute_properties,
            'data': compute_data,
        }

        signal(SIGINT, SIG_DFL)

        self.mpi_manager = MPIManager()

        # initialize mol
        self.mol = Molecule(project, input_file, log, silent=silent)
        self.mol.usempi = usempi

        if input_dict:
            self.mol.load_config(input_dict)
        else:
            self.mol.load_config(input_file)
        if self.mpi_manager.rank != 0:
            if os.name == 'nt':  # Windows
                log = 'NUL'
            else:
                log = '/dev/null'
        else:
            log = self.mol.log

        self.mol.data["OQP::log_filename"] = log

        # Apply the parsed `omp_threads` at runtime. The pre-import hook only sees
        # a CLI flag or an input *file*, so this also covers the programmatic
        # Runner(input_dict=...) path. omp_set_num_threads takes effect for the
        # subsequent SCF parallel regions; we also export OMP_NUM_THREADS so the
        # input checker / any later BLAS sizing see a consistent value.
        _omp = self.mol.config.get("input", {}).get("omp_threads", 0)
        if _omp and _omp > 0:
            if oqp.lib.oqp_have_openmp():
                oqp.lib.oqp_omp_set_num_threads(int(_omp))
                os.environ["OMP_NUM_THREADS"] = str(int(_omp))
            elif not self.mol.silent:
                print(f"PyOQP WARNING: omp_threads={_omp} requested but this "
                      "OpenQP build has no OpenMP support; running serially.")
        # check input values set default omp_num_threads
        _input_file = getattr(self.mol, "input_file", None)
        check_input_values(
            self.mol.config,
            input_dir=os.path.dirname(os.path.abspath(_input_file)) if _input_file else None,
        )

        # Reload mol if possible
        self.reload()

        # Attach the starting time to mol
        self.mol.start_time = start_time
        # Set up banner
        oqp.oqp_banner(self.mol)

        dump_log(self.mol, title='', section='start')
        dump_log(self.mol, title='PyOQP: Symmetry metadata', section='symmetry')
        # Log the running OpenQP version + git HEAD commit so every output records
        # which build produced it (works for source/editable installs; falls back to
        # the package version when no git tree is reachable).
        try:
            import subprocess as _sp
            _pkg = os.path.dirname(os.path.abspath(oqp.__file__))
            # This very module: present in any OpenQP package, used as a sentinel.
            _sentinel = os.path.join(_pkg, "pyoqp.py")
            _ver = getattr(oqp, "__version__", "")
            _commit = ""
            for _d in (os.environ.get("OPENQP_ROOT"), _pkg):
                if not _d:
                    continue
                # Only trust a discovered git worktree if it actually TRACKS the
                # OpenQP package source. A wheel installed under an unrelated
                # project's virtualenv would otherwise let `git -C <site-packages>/oqp`
                # discover that parent project's .git and report its commit as the
                # OpenQP build. Requiring this package file to be tracked rules that
                # out -- it is tracked only in a genuine source/editable checkout.
                _chk = _sp.run(["git", "-C", _d, "ls-files", "--error-unmatch", _sentinel],
                               capture_output=True, text=True, timeout=2)
                if _chk.returncode != 0:
                    continue
                _r = _sp.run(["git", "-C", _d, "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=2)
                if _r.returncode == 0 and _r.stdout.strip():
                    _commit = _r.stdout.strip()
                    # "+dirty" reflects only the OpenQP package subtree.
                    _dirty = _sp.run(["git", "-C", _d, "status", "--porcelain", "--", _pkg],
                                     capture_output=True, text=True, timeout=2)
                    if _dirty.stdout.strip():
                        _commit += "+dirty"
                    break
            _bv = "OpenQP" + (" v%s" % _ver if _ver else "") + \
                  ((" (git %s)" % _commit) if _commit else " (git commit unknown)")
            dump_log(self.mol, section='', title='PyOQP build: %s' % _bv)
        except Exception:
            pass

    def run(self, test_mod=False):
        """
        Run the OQP calculation or test.

        Args:
            test_mod (bool): Flag to run in test mode.
        """

        # Get the run type from mol configuration
        run_type = self.mol.config["input"]["runtype"]

        # Turn off convergence check for test
        if test_mod:
            self.mol.config['tests']['exception'] = True

        # Run calculations
        self.run_func[run_type](self.mol)

        dump_log(self.mol, title='', section='end')

    def results(self):
        """ Collect calculation results for internal Python communication.

        Returns:
            dict: A dictionary containing calculation results.
        """
        summary = {
            'atoms': self.mol.get_atoms(),
            'system': self.mol.get_system(),
            'energy': self.mol.energies,
            'grad': self.mol.grads,
            'dcm': self.mol.dcm,
            'nac': self.mol.nac,
            'soc': self.mol.soc,
            'data': self.mol.get_data(),
        }
        return summary

    def reload(self):
        """
        Reload calculation based on the specified guess type.
        """
        guess_type = self.mol.config['guess']['type']
        guess_file = self.mol.config['guess']['file']

        if guess_type == 'json':
            self.mol.load_data()
        elif guess_type == 'auto':
            if os.path.exists(guess_file):
                self.mol.load_data()
        else:
            pass

    def back_door(self, data):
        """
        Set back door data for the molecule.

        Args:
            data: Data to be set as back door for the molecule.
        """
        self.mol.back_door = data

    def test(self):
        """
        Run tests and check results against reference data (json).
        """
        if self.mpi_manager.rank == 0:
            message, diff = self.mol.check_ref()
        else:
            message = None
            diff = None
        return message, diff

def _warn_if_no_openmp():
    """Warn when threads were requested but liboqp was built without OpenMP, so
    the request (input omp_threads / --omp / OMP_NUM_THREADS) has no effect."""
    try:
        want = int(os.environ.get("OMP_NUM_THREADS", "1"))
    except ValueError:
        want = 1
    if want <= 1:
        return
    try:
        have = bool(oqp.lib.oqp_have_openmp())
    except Exception:
        return
    if not have:
        print(f"PyOQP WARNING: {want} OpenMP threads were requested "
              "(omp_threads/--omp/OMP_NUM_THREADS) but this OpenQP build has no "
              "OpenMP support; the calculation will run serially.")


def main():
    """
    Main function to handle command-line arguments and run OQP.
    """
    _warn_if_no_openmp()
    parser = argparse.ArgumentParser(description='OQP Runner',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input', nargs='?', help='Input file')
    parser.add_argument('--run_tests', '--test',
                        dest='run_tests',
                        metavar='path',
                        help='run tests from a specified folder or:\n'
                             '  all    - Run all tests in examples\n'
                             '  other  - Run tests in examples/other')
    parser.add_argument('--silent', action='store_true', help='run silently')
    parser.add_argument('--nompi', action='store_true', help='disable mpi functions')
    parser.add_argument('--omp', metavar='N', type=int,
                        help='OpenMP threads per process/MPI rank (overrides the\n'
                             "input's omp_threads and OMP_NUM_THREADS; applied\n"
                             'before the OpenMP runtime loads)')
    args = parser.parse_args()

    if args.run_tests:
        report, status = run_tests(args.run_tests)
        print(report)
        sys.exit(status)
    if not args.input:
        parser.print_help()
        sys.exit(1)

    input_file = os.path.abspath(args.input)
    if not os.path.exists(input_file):
        print(f'\n   PyOQP input file {input_file} not found\n')
        sys.exit(1)

    input_path = os.path.dirname(input_file)
    project_name = os.path.basename(input_file).replace('.inp', '')
    log = f'{input_path}/{project_name}.log'

    mpi_manager = MPIManager()
    usempi = True if not args.nompi and mpi_manager.use_mpi > 0 else False

    if usempi:
        input_file = mpi_manager.bcast(input_file)
        project_name = mpi_manager.bcast(project_name)
        log = mpi_manager.bcast(log)

    silent = 1 if args.silent else 0

    # Detect the OpenMM-based QM/MM-MD mode without importing OpenMM-dependent
    # modules (so plain energy/grad/NAMD runs work without OpenMM installed).
    qmmm_flag = False
    try:
        import configparser
        _cfg = configparser.ConfigParser()
        _cfg.read(input_file)
        qmmm_flag = _cfg.getboolean('input', 'qmmm_flag', fallback=False)
    except Exception:
        qmmm_flag = False
    # runtype=namd with qmmm goes through the Runner (-> compute_namd ->
    # NAMD_QMMM); only the ground-state OpenMM-MD path is handled here.
    runtype_l = ""
    try:
        runtype_l = _cfg.get('input', 'runtype', fallback='energy').strip().lower()
    except Exception:
        runtype_l = 'energy'
    if qmmm_flag and runtype_l != 'namd':
        from oqp.library.qmmm_md import QMMM_MD
        md = QMMM_MD(oqp_cfg=input_file)
        md.run()
        return

    # Initialize OQP class
    oqp_runner = Runner(project=project_name,
                        input_file=input_file,
                        log=log,
                        silent=silent,
                        usempi=usempi,
                        )

    # Run OQP
    oqp_runner.run()
    oqp_runner.results()
    mpi_manager.finalize_mpi()


def run_tests(test_path):
    """
    Run OQP tests.

    Args:
        test_path (str): Path to the test directory or 'all' or 'other'.

    Returns:
        str: Test report.
    """
    from oqp.utils.oqp_tester import OQPTester
    mpi_manager = MPIManager()
    tester = OQPTester(base_test_dir=None,
                       output_dir='openqp_test_tmp',
                       total_cpus=None,
                       omp_threads=4, mpi_manager=mpi_manager)
    return tester.run(test_path), tester.status


if __name__ == "__main__":
    main()
