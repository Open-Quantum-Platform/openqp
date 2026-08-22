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


def _strip_preimport_comments(text):
    """Strip .inp/.oqp comments without importing the semantic parser."""
    cleaned = []
    for line in text.splitlines():
        quote = None
        escaped = False
        kept = []
        for index, char in enumerate(line):
            if escaped:
                kept.append(char)
                escaped = False
                continue
            if char == "\\" and quote:
                kept.append(char)
                escaped = True
                continue
            if char in {'"', "'"}:
                quote = None if quote == char else (char if quote is None else quote)
                kept.append(char)
                continue
            if quote is None and char in {"#", "!"} and (
                    index == 0 or line[index - 1].isspace()):
                break
            kept.append(char)
        cleaned.append("".join(kept))
    return "\n".join(cleaned)


def _apply_omp_threads_from_input(argv):
    """Honour an OpenMP thread-count request from the input file or CLI BEFORE
    the native OpenMP runtime initialises (it caches OMP_NUM_THREADS when liboqp
    loads, so this must run before `import oqp`).

    Sources, highest precedence first:
      * CLI:   --omp N   (or --omp=N)
      * input: a line  `omp_threads = N`  (typically in the [input] section)

    The value sets OMP_NUM_THREADS (threads per process / MPI rank).  For a
    native DFTB input, it also seeds MKL_NUM_THREADS when the user has not set
    it explicitly: the DFTB eigensolver manages its dense-MKL region locally,
    while its Davidson sigma columns remain pinned to one BLAS thread.  If
    neither source is given, the existing environment / built-in default is
    left untouched.
    """
    import re
    n = None
    input_text = ""
    for i, a in enumerate(argv):
        if a in ("--omp", "--omp-threads") and i + 1 < len(argv):
            n = argv[i + 1]
            break
        if a.startswith("--omp="):
            n = a.split("=", 1)[1]
            break
    inp = next((a for a in argv[1:]
                if not a.startswith("-") and os.path.isfile(a)), None)
    if inp:
        try:
            with open(inp, encoding="utf-8", errors="ignore") as fh:
                input_text = _strip_preimport_comments(fh.read())
            if n is None:
                # Legacy INI, canonical .oqp section/top-level options, and a
                # small controlled-natural-language spelling.  This scan must
                # stay dependency-free because it runs before liboqp/OpenMP is
                # imported; the full .oqp parser validates the value later.
                patterns = (
                    r"(?mi)^[ \t]*omp_threads[ \t]*=[ \t]*(\d+)",
                    r"(?i)\bomp_threads[ \t]*=[ \t]*(\d+)",
                    r"(?i)\b(?:use|with|using)?[ \t]*(\d+)[ \t]+(?:openmp[ \t]+)?threads?\b",
                    r"(?i)(\d+)[ \t]*(?:개[ \t]*)?(?:스레드|쓰레드)",
                )
                for pattern in patterns:
                    m = re.search(pattern, input_text)
                    if m:
                        n = m.group(1)
                        break
        except OSError:
            pass
    if n is not None:
        try:
            ni = int(n)
        except (TypeError, ValueError):
            return
        if ni >= 1:
            os.environ["OMP_NUM_THREADS"] = str(ni)
            native_dftb = bool(re.search(r"(?i)\b[a-z0-9-]*dftb\b", input_text))
            if native_dftb:
                # Respect an explicit BLAS policy.  When it is absent, seed
                # MKL before the runtime loads so the dense reference
                # eigensolver can use the requested cores.  The DFTB library
                # locally serializes narrow Davidson sigma columns, avoiding
                # nested OMP x MKL oversubscription in the response stage.
                os.environ.setdefault("MKL_NUM_THREADS", str(ni))


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


# Honour explicit per-rank threading before populating conservative defaults:
# this lets an input request seed the DFTB/MKL dense-eigensolver policy while
# preserving every user-exported thread setting.
_apply_omp_threads_from_input(sys.argv)
_set_threading_defaults()

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


def _resolve_semantic_input(input_file, mpi_manager=None, usempi=None):
    """Resolve a ``.oqp`` file once and return an MPI-safe execution payload.

    Only rank zero may write a correction ``*.resolved.oqp`` artifact.  Prose
    is never executed; canonical requests are lowered and broadcast so every
    rank runs exactly the same validated configuration.
    """
    from oqp.utils.oqp_input import resolve_oqp_file

    mpi_manager = mpi_manager or MPIManager()
    broadcast = bool(mpi_manager.use_mpi if usempi is None else
                     (usempi and mpi_manager.use_mpi))
    envelope = None
    if mpi_manager.rank == 0 or not broadcast:
        try:
            resolution = resolve_oqp_file(input_file, write_resolved=True)
            source = str(resolution.source_path or os.path.abspath(input_file))
            effective = str(resolution.resolved_path or resolution.source_path or input_file)
            payload = {
                "source": source,
                "resolved": effective,
                "was_natural": bool(resolution.was_natural),
                "canonical_text": resolution.canonical_text,
                "legacy_config": {
                    section: dict(values)
                    for section, values in resolution.legacy_config.items()
                },
            }
            envelope = {"ok": True, "payload": payload}
        except Exception as err:
            envelope = {"ok": False, "error": str(err)}
    if broadcast:
        envelope = mpi_manager.bcast(envelope)
    if not envelope["ok"]:
        from oqp.utils.oqp_input import OQPInputError
        raise OQPInputError(envelope["error"])
    return envelope["payload"]


def _flatten_config(config):
    """Flatten a section dictionary for the legacy ground-state QM/MM driver."""
    return {
        "%s.%s" % (section, key): value
        for section, values in config.items()
        for key, value in values.items()
    }


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


# The OMP_NUM_THREADS in force before any Runner touched it.  Captured once so a
# Runner with an explicit [input] omp_threads cannot change what a later default
# Runner sees.
_RUNNER_OMP_BASELINE = [os.environ.get("OMP_NUM_THREADS")]


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
                 input_dict=None, silent=0, usempi=True, input_metadata=None):
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

        self.mpi_manager = MPIManager()

        # Programmatic Runner(input_file="request.oqp") follows the same path as
        # the CLI.  Supplying input_dict explicitly retains its historical
        # precedence and does not reinterpret the named file.
        if input_dict is None and input_file and str(input_file).lower().endswith('.oqp'):
            input_metadata = _resolve_semantic_input(
                input_file, self.mpi_manager, usempi=usempi)
            if input_metadata["was_natural"]:
                raise ValueError(
                    "Free-form text is not a runnable .oqp input. Review the corrected "
                    "canonical file at %s and run that file instead."
                    % input_metadata["resolved"]
                )
            input_file = input_metadata["resolved"]
            input_dict = input_metadata["legacy_config"]

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

        # initialize mol
        self.mol = Molecule(project, input_file, log, silent=silent)
        self.mol.usempi = usempi
        # NAMD keeps backward compatibility with legacy [md] triplet labels
        # (T1 is the first root), while canonical .oqp files use public labels
        # (T0 is the first root).  Preserve that provenance after lowering.
        self.mol.oqp_public_state_labels = bool(input_metadata)
        if input_metadata:
            self.mol.oqp_input_source = input_metadata.get("source")
            self.mol.oqp_resolved_input = input_metadata.get("resolved")
            self.mol.oqp_input_was_natural = bool(input_metadata.get("was_natural"))
            self.mol.oqp_canonical_input = input_metadata.get("canonical_text", "")

        if input_dict is not None:
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

        # (The `perf` preset is resolved inside mol.load_config, before the config is
        #  pushed to the control struct; mol.perf_report holds the resolution for logging.)

        # Apply the parsed `omp_threads` at runtime. The pre-import hook only sees
        # a CLI flag or an input *file*, so this also covers the programmatic
        # Runner(input_dict=...) path. omp_set_num_threads takes effect for the
        # subsequent SCF parallel regions; we also export OMP_NUM_THREADS so the
        # input checker / any later BLAS sizing see a consistent value.
        _omp = self.mol.config.get("input", {}).get("omp_threads", 0)
        if not (_omp and _omp > 0):
            # A previous Runner in this process may have set an explicit count,
            # the environment variable and OMP_THREADS_FROM_ENV.  None of that
            # was ever undone, so a later default job inherited the earlier
            # job's thread count and _fci_lib_threads treated the stale value as
            # an explicit request.  Restore the default policy for this Runner.
            oqp.OMP_THREADS_FROM_ENV = False
            if getattr(self, "_omp_default_env", None) is None:
                type(self)._omp_default_env = os.environ.get("OMP_NUM_THREADS")
            if _RUNNER_OMP_BASELINE[0] is not None:
                os.environ["OMP_NUM_THREADS"] = _RUNNER_OMP_BASELINE[0]
                if oqp.lib.oqp_have_openmp():
                    oqp.lib.oqp_omp_set_num_threads(int(_RUNNER_OMP_BASELINE[0]))
        if _omp and _omp > 0:
            if oqp.lib.oqp_have_openmp():
                oqp.lib.oqp_omp_set_num_threads(int(_omp))
                os.environ["OMP_NUM_THREADS"] = str(int(_omp))
                # The input asked for this count, so the kernels that take an
                # explicit thread argument must use it rather than their own
                # default (see oqp.library.fci._fci_lib_threads).
                oqp.OMP_THREADS_FROM_ENV = True
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
        # Initialize the log before the native banner appends to it.  The old
        # order wrote the banner first and immediately truncated it in the
        # ``start`` section, hiding the contributor and resource information.
        dump_log(self.mol, title='', section='start',
                 info={"build": _openqp_build_label()})
        oqp.oqp_banner(self.mol)
        dump_log(self.mol, title='PyOQP: Calculation request', section='calculation')
        dump_log(self.mol, title='PyOQP: Symmetry metadata', section='symmetry')
        self._log_perf_settings()

    def _log_perf_settings(self):
        """Append the resolved performance settings + warnings to the log."""
        report = getattr(self.mol, "perf_report", None)
        if not report:
            return
        from oqp.utils import perf_levels
        block = perf_levels.format_report(getattr(self.mol, "perf_level", perf_levels.UNSET),
                                          report, getattr(self.mol, "perf_warns", []))
        if not block:
            return
        if getattr(self.mol, "log", None):
            try:
                with open(self.mol.log, 'a', encoding='utf-8') as fout:
                    fout.write(block + "\n")
            except OSError:
                pass
        if not getattr(self.mol, "silent", False):
            print(block)

    def run(self, test_mod=False):
        """
        Run the OQP calculation or test.

        Args:
            test_mod (bool): Flag to run in test mode.
        """

        # A reused Runner/Molecule carries the previous calculation's gradient
        # buffer AND its validity flag.  The flag is only initialised in
        # Molecule.__init__, so a gradient run followed by an energy-only run on
        # the same object would publish the earlier gradient as if this run had
        # produced it -- the stale-result failure #336 is about, but with a
        # physically plausible value instead of uninitialised memory.  Any
        # multi-step consumer that reuses one Molecule (OPENQP.set/run,
        # optimisation, NAMD) is exposed.  Invalidate at the start of every
        # top-level calculation; the kernels re-declare it when they write.
        self.mol._grad_valid = False

        # Get the run type from mol configuration
        run_type = self.mol.config["input"]["runtype"]

        # Turn off convergence check for test
        if test_mod:
            self.mol.config['tests']['exception'] = True

        qmmm_flag = self.mol.config["input"].get("qmmm_flag", False)
        qmmm_flag = (qmmm_flag is True or
                     str(qmmm_flag).strip().lower() in {"true", "1", "yes", "on"})

        # The OpenMM QM/MM-MD driver is a distinct workflow from the all-QM
        # velocity-Verlet driver.  The CLI already dispatches it before Runner;
        # keep the advertised programmatic Runner(input_file="job.oqp") path
        # equivalent, and cover legacy programmatic inputs at the same time.
        if qmmm_flag and str(run_type).strip().lower() == "md":
            if self.mol.usempi and self.mpi_manager.use_mpi > 0:
                raise RuntimeError(
                    "ground-state QM/MM-MD cannot run under MPI; create Runner with usempi=False"
                )
            from oqp.library.qmmm_md import QMMM_MD
            self.qmmm_md = QMMM_MD(mol=self.mol)
            self.qmmm_md.run()
        else:
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
        state_tracking = self.mol.get_state_tracking()
        if state_tracking is not None:
            summary['state_tracking'] = state_tracking
        return summary

    def reload(self):
        """
        Reload calculation based on the specified guess type.
        """
        input_config = self.mol.config.get('input', {})
        md_config = self.mol.config.get('md', {})
        restart_namd = (
            str(input_config.get('runtype', '')).strip().lower() == 'namd'
            and (
                md_config.get('restart') is True
                or str(md_config.get('restart', '')).strip().lower()
                in {'true', '1', 'yes', 'on'}
            )
        )
        if restart_namd:
            # The atomic NAMD checkpoint is the authoritative electronic
            # state. In particular, a save_mol JSON file is a mutable output
            # and may be absent or incomplete after interruption; preloading
            # it would fail before NAMD can validate and restore the checkpoint.
            return
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
    parser.add_argument('input', nargs='?', help='Input file (.inp or .oqp)')
    parser.add_argument('--run_tests', '--run-tests', '--test',
                        dest='run_tests',
                        metavar='path',
                        help='run tests from a specified folder or:\n'
                             '  all    - Run the standard suite in examples\n'
                             '           (slow/non-self-contained exclusions apply)\n'
                             '  other  - Run tests in examples/other\n'
                             '  WF_methods - Run curated wavefunction examples')
    parser.add_argument('--input-format', '--input_format', '--test-inputs',
                        dest='test_input_format',
                        choices=('auto', 'inp', 'oqp', 'both'),
                        help='input syntax selected by --run_tests:\n'
                             '  auto   - prefer .oqp; keep representative legacy .inp (default)\n'
                             '  inp    - select legacy .inp\n'
                             '  oqp    - select concise .oqp\n'
                             '  both   - include .inp and .oqp in the selected test scope')
    parser.add_argument('--validate_examples', dest='validate_examples',
                        metavar='dir', nargs='?', const='',
                        help='validate that every example reference under DIR\n'
                             '(default: $OPENQP_ROOT/share/examples) carries the\n'
                             'regression values its runtype requires, then exit')
    parser.add_argument('--check_feature_coverage', dest='check_feature_coverage',
                        metavar='dir', nargs='?', const='',
                        help='check that every opt-in feature flag in the input\n'
                             'schema is exercised by at least one example under DIR\n'
                             '(default: $OPENQP_ROOT/share/examples); fails if a new\n'
                             'feature ships without a test')
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

    if args.test_input_format is not None and not args.run_tests:
        parser.error('--input-format is only valid with --run_tests')

    if args.generate_reference:
        sys.exit(generate_reference_cli(args.generate_reference))

    if args.check_feature_coverage is not None:
        sys.exit(feature_coverage_cli(args.check_feature_coverage))

    if args.validate_examples is not None:
        sys.exit(validate_examples_cli(args.validate_examples))

    if args.run_tests:
        try:
            report, status = run_tests(
                args.run_tests,
                input_format=args.test_input_format or 'auto',
            )
        except ValueError as err:
            parser.error(str(err))
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
    project_name = os.path.splitext(os.path.basename(input_file))[0]
    log = f'{input_path}/{project_name}.log'

    mpi_manager = MPIManager()
    usempi = True if not args.nompi and mpi_manager.use_mpi > 0 else False

    if usempi:
        input_file = mpi_manager.bcast(input_file)
        project_name = mpi_manager.bcast(project_name)
        log = mpi_manager.bcast(log)

    silent = 1 if args.silent else 0

    semantic_input = None
    input_dict = None
    if input_file.lower().endswith('.oqp'):
        try:
            semantic_input = _resolve_semantic_input(
                input_file, mpi_manager, usempi=usempi)
        except (OSError, ValueError) as err:
            if mpi_manager.rank == 0:
                print(f'\n   PyOQP .oqp input error: {err}\n')
            if usempi:
                mpi_manager.finalize_mpi()
            sys.exit(2)
        if semantic_input["was_natural"]:
            if mpi_manager.rank == 0:
                print('\n   PyOQP: free-form text was not executed.\n'
                      '   A corrected canonical .oqp suggestion was written to:\n'
                      f'   {semantic_input["resolved"]}\n'
                      '   Review it, then run that canonical file.\n')
            if usempi:
                mpi_manager.finalize_mpi()
            sys.exit(2)
        input_file = semantic_input["resolved"]
        input_dict = semantic_input["legacy_config"]

    # Detect the OpenMM-based QM/MM-MD mode without importing OpenMM-dependent
    # modules (so plain energy/grad/NAMD runs work without OpenMM installed).
    qmmm_flag = False
    _cfg = None
    if input_dict is not None:
        qmmm_flag = str(input_dict.get('input', {}).get('qmmm_flag', 'False')).lower() \
            in {'true', '1', 'yes', 'on'}
    else:
        try:
            import configparser
            _cfg = configparser.ConfigParser()
            _cfg.read(input_file)
            qmmm_flag = _cfg.getboolean('input', 'qmmm_flag', fallback=False)
        except Exception:
            qmmm_flag = False
    # Only the legacy ground-state OpenMM-MD path is handled here. QM/MM energy,
    # gradient, optimization, and NAMD inputs go through Runner so their normal
    # runtype dispatch and read_system paths are preserved.
    runtype_l = ""
    if input_dict is not None:
        runtype_l = str(input_dict.get('input', {}).get('runtype', 'energy')).strip().lower()
    else:
        try:
            runtype_l = _cfg.get('input', 'runtype', fallback='energy').strip().lower()
        except Exception:
            runtype_l = 'energy'
    if qmmm_flag and runtype_l == 'md':
        if usempi:
            if mpi_manager.rank == 0:
                print('\n   PyOQP: ground-state QM/MM-MD cannot run under MPI; '
                      'rerun with --nompi.\n')
            mpi_manager.finalize_mpi()
            sys.exit(2)
        from oqp.library.qmmm_md import QMMM_MD
        md = QMMM_MD(oqp_cfg=_flatten_config(input_dict) if input_dict is not None
                     else input_file)
        md.run()
        return

    # Initialize OQP class
    oqp_runner = Runner(project=project_name,
                        input_file=input_file,
                        log=log,
                        input_dict=input_dict,
                        input_metadata=semantic_input,
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
    identity/metadata) next to the input, dropping internal OQP:: arrays. This is
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
        ref = os.path.splitext(inp)[0] + '.json'
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
        ref = os.path.splitext(inp)[0] + '.json'
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


def _default_examples_dir():
    """$OPENQP_ROOT/share/examples (the installed runtime examples)."""
    root = os.environ.get('OPENQP_ROOT') or getattr(
        oqp, 'oqp_root', os.path.dirname(os.path.abspath(__file__)))
    examples_dir = os.path.join(root, 'share', 'examples')
    if not os.path.isdir(examples_dir):
        examples_dir = os.path.join(root, 'examples')
    return examples_dir


def feature_coverage_cli(examples_dir):
    """Gate: every opt-in feature flag in the input schema must be exercised by
    at least one example, so a new feature cannot ship without a test. Exit 0
    if covered (grandfathered gaps are warnings), 1 if a new flag is untested."""
    from oqp.utils import regression
    examples_dir = examples_dir or _default_examples_dir()
    failures, grandfathered = regression.feature_coverage(examples_dir)
    for flag, why in grandfathered:
        print(f'   PyOQP NOTE: feature {flag} has no example yet (tracked gap): {why}')
    if not failures:
        print(f'PyOQP: every opt-in feature flag under {examples_dir} is '
              f'exercised by an example (or classified).')
        return 0
    print(f'PyOQP: {len(failures)} feature flag(s) ship without a test '
          f'(add an example, or classify in oqp/utils/regression.py):')
    for flag, why in failures:
        print(f'   {flag:<34} {why}')
    return 1


def validate_examples_cli(examples_dir):
    """Gate: every example reference must carry the regression values its
    runtype/method/properties require (per the registry). Returns a process
    exit code (0 = all good, 1 = at least one reference is missing a value)."""
    from oqp.utils import regression
    if not examples_dir:
        # Default to the installed runtime root's examples -- the same location
        # OQPTester uses (oqp.oqp_root/share/examples). NB: pyoqp.py lives inside
        # the oqp package, so do NOT go up a level (that lands on site-packages
        # and the gate then scans nothing).
        root = os.environ.get('OPENQP_ROOT') or getattr(
            oqp, 'oqp_root', os.path.dirname(os.path.abspath(__file__)))
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


def run_tests(test_path, *, input_format='auto'):
    """
    Run OQP tests.

    Args:
        test_path (str): Path to the test directory or 'all' or 'other'.
        input_format (str): ``auto`` (default), ``inp``, ``oqp``, or ``both``.

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
    return tester.run(test_path, input_format=input_format), tester.status


if __name__ == "__main__":
    main()
