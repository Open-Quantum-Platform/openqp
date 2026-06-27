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
        "OMP_STACKSIZE": "256M",
        "GOMP_STACKSIZE": "256M",
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
   compute_energy, compute_grad, compute_nac, compute_soc, compute_geom,
   compute_nacme, compute_properties, compute_data, compute_hess, compute_thermo
)
from oqp.utils.mpi_utils import MPIManager


def _openqp_build_label():
    """Return a concise OpenQP package/version label for the first log block."""
    try:
        import subprocess as _sp
        pkg_dir = os.path.dirname(os.path.abspath(oqp.__file__))
        sentinel = os.path.join(pkg_dir, "pyoqp.py")
        version = getattr(oqp, "__version__", "")
        commit = ""

        for candidate in (getattr(oqp, "oqp_root", None), pkg_dir):
            if not candidate:
                continue
            root = _sp.run(["git", "-C", candidate, "rev-parse", "--show-toplevel"],
                           capture_output=True, text=True, timeout=2)
            if root.returncode != 0 or not root.stdout.strip():
                continue
            root_dir = root.stdout.strip()
            sentinel_rel = os.path.relpath(sentinel, root_dir)
            tracked = _sp.run(["git", "-C", root_dir, "ls-files", "--error-unmatch", sentinel_rel],
                              capture_output=True, text=True, timeout=2)
            if tracked.returncode != 0:
                continue
            head = _sp.run(["git", "-C", root_dir, "rev-parse", "HEAD"],
                           capture_output=True, text=True, timeout=2)
            if head.returncode == 0 and head.stdout.strip():
                commit = head.stdout.strip()
                pkg_rel = os.path.relpath(pkg_dir, root_dir)
                dirty = _sp.run(["git", "-C", root_dir, "status", "--porcelain", "--", pkg_rel],
                                capture_output=True, text=True, timeout=2)
                if dirty.stdout.strip():
                    commit += "+dirty"
                break

        return "OpenQP" + (" v%s" % version if version else "") + \
            ((" (git HEAD %s)" % commit) if commit else " (git HEAD unknown)")
    except Exception:
        version = getattr(oqp, "__version__", "")
        return "OpenQP" + (" v%s" % version if version else "") + " (git HEAD unknown)"


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

        dump_log(self.mol, title='', section='start',
                 info={"build": _openqp_build_label()})
        dump_log(self.mol, title='PyOQP: Symmetry metadata', section='symmetry')

    def run(self, test_mod=False):
        """
        Run the OQP calculation or test.

        Args:
            test_mod (bool): Flag to run in test mode.
        """
        # Set up logfile

        if self.mpi_manager.rank != 0:
            if os.name == 'nt':  # Windows
                log = 'NUL'
            else:
                log = '/dev/null'
        else:
            log = self.mol.log

        self.mol.data["OQP::log_filename"] = log

        # Set up banner
        oqp.oqp_banner(self.mol)

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
    parser.add_argument('--validate_examples', dest='validate_examples',
                        metavar='dir', nargs='?', const='',
                        help='validate that every example reference under DIR\n'
                             '(default: $OPENQP_ROOT/share/examples) carries the\n'
                             'regression values its runtype requires, then exit')
    parser.add_argument('--generate_reference', dest='generate_reference',
                        metavar='input.inp', nargs='+',
                        help='run each INPUT and (re)write its lean .json test\n'
                             'reference next to it -- only the regression-registry\n'
                             'keys are kept (internal OQP:: arrays dropped). This\n'
                             'is the supported way to add/refresh example references.')
    parser.add_argument('--silent', action='store_true', help='run silently')
    parser.add_argument('--nompi', action='store_true', help='disable mpi functions')
    parser.add_argument('--omp', metavar='N', type=int,
                        help='OpenMP threads per process/MPI rank (overrides the\n'
                             "input's omp_threads and OMP_NUM_THREADS; applied\n"
                             'before the OpenMP runtime loads)')
    args = parser.parse_args()

    if args.generate_reference:
        sys.exit(generate_reference_cli(args.generate_reference))

    if args.validate_examples is not None:
        sys.exit(validate_examples_cli(args.validate_examples))

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


def generate_reference_cli(inputs):
    """(Re)generate lean JSON test references for the given input files.

    Runs each input and writes only the regression-registry keys (physics +
    identity/metadata) next to the .inp, dropping internal OQP:: arrays. This is
    the supported, repeatable way to add or refresh example references so the
    committed set stays clean, small, and consistent with the registry. The
    written references are validated against the registry before returning.
    """
    import tempfile
    from oqp.utils import regression

    rc = 0
    for inp in inputs:
        inp = os.path.abspath(inp)
        if not os.path.exists(inp):
            print(f'   PyOQP generate_reference: input {inp} not found')
            rc = 1
            continue
        project = os.path.splitext(os.path.basename(inp))[0]
        ref = inp[:-4] + '.json' if inp.endswith('.inp') else inp + '.json'
        with tempfile.TemporaryDirectory() as tmp:
            runner = Runner(project=project, input_file=inp,
                            log=os.path.join(tmp, project + '.log'),
                            silent=1, usempi=False)
            runner.run(test_mod=True)
            # Point save_data at the example location and write the lean bundle.
            runner.mol.log = ref.replace('.json', '.log')
            runner.mol.save_data(lean=True)
        print(f'   PyOQP wrote lean reference {ref}')

    # Validate what we just wrote.
    failures = []
    for inp in inputs:
        inp = os.path.abspath(inp)
        ref = inp[:-4] + '.json' if inp.endswith('.inp') else inp + '.json'
        if not os.path.exists(ref):
            continue
        runtype, excited, props = regression._context_from_input(inp)
        miss = regression.missing_required(
            regression._present_nonempty(ref, inp), runtype, excited, props)
        if miss:
            failures.append((inp, miss))
    for inp, miss in failures:
        print(f'   PyOQP WARNING: {inp} reference still missing: {", ".join(miss)}')
        rc = 1
    return rc


def validate_examples_cli(examples_dir):
    """Gate: every example reference must carry the regression values its
    runtype/method/properties require (per the registry). Returns a process
    exit code (0 = all good, 1 = at least one reference is missing a value)."""
    from oqp.utils import regression
    if not examples_dir:
        root = os.environ.get('OPENQP_ROOT') or os.path.abspath(
            os.path.join(os.path.dirname(__file__), os.pardir))
        examples_dir = os.path.join(root, 'share', 'examples')
        if not os.path.isdir(examples_dir):
            examples_dir = os.path.join(root, 'examples')
    failures = regression.validate_examples(examples_dir)
    if not failures:
        print(f'PyOQP: all example references under {examples_dir} carry their '
              f'required regression values.')
        return 0
    print(f'PyOQP: {len(failures)} example reference(s) missing required '
          f'regression values (see the registry in oqp/utils/regression.py):')
    for inp, miss in failures:
        print(f'   {inp:<60} missing: {", ".join(miss)}')
    return 1


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
    # OMP threads per test. Fewer threads -> more tests run concurrently
    # (max_workers = total_cpus // omp_threads), which is much faster on the
    # small CI examples; override with OQP_TEST_OMP_THREADS.
    try:
        omp_threads = max(1, int(os.environ.get("OQP_TEST_OMP_THREADS", "4")))
    except ValueError:
        omp_threads = 4
    tester = OQPTester(base_test_dir=None,
                       output_dir='openqp_test_tmp',
                       total_cpus=None,
                       omp_threads=omp_threads, mpi_manager=mpi_manager)
    return tester.run(test_path), tester.status


if __name__ == "__main__":
    main()
