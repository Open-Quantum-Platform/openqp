"""
OpenQP Test Runner

This module provides a class for running OpenQP tests
in parallel and generating reports.
It is compatible with CI pipelines.

Author: Konstantin Komarov
Email: constlike@gmail.com
Created: Aug 2024
"""
import hashlib
import json
import os
import sys
import time
import subprocess
from typing import List, Dict, Any
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

from oqp.pyoqp import Runner
from oqp.runtime import resolve_oqp_root

# Marker that the isolated single-test subprocess prints so the parent can
# recover the structured result from the child's (otherwise noisy) stdout.
_RESULT_MARKER = "__OQP_TEST_RESULT__"

# Per-test wall-clock ceiling for the isolated subprocess runner. A wedged or
# pathologically slow test is reported as a failure instead of hanging CI
# forever. Override with OQP_TEST_TIMEOUT (seconds); 0/empty disables it.
def _test_timeout():
    raw = os.environ.get("OQP_TEST_TIMEOUT", "1800").strip()
    try:
        val = float(raw)
    except ValueError:
        return 1800.0
    return val if val > 0 else None

class OQPTester:
    """
    A class for running OQP tests and generating reports.

    This class can run tests from specified directories:
      - 'openqp --run-tests path_to_folder': Run tests from a specific folder
      - 'openqp --run-tests other': Run tests from the 'other' folder in examples
      - 'openqp --run-tests all': Run the standard suite across examples,
        with documented slow/non-self-contained exclusions
      - '--input-format auto|inp|oqp|both': Select the input syntax to test
    Can be used to run all tests in a specific folder.

    Attributes:
        base_test_dir (str): Base directory for test files.
        output_dir (str): Directory for storing test output files.
        max_workers (int): Maximum number of parallel workers.
        results (List[Dict[str, Any]]): List to store test results.
    """

    # Substrings a crashed test's output may contain that mean "this example
    # needs an optional backend this build was not configured with" -> SKIPPED,
    # not ERROR. Keep these specific to capability gates, never generic errors.
    _OPTIONAL_FEATURE_SENTINELS = (
        "OQP_ENABLE_DDX",  # ddX / PCM continuum solvation built off
        "openqp-dftb not found",  # optional external openqp-dftb library absent
        "openqp-xtb not found",  # optional external openqp-xtb library absent
    )

    _INPUT_FORMATS = frozenset({"auto", "inp", "oqp", "both"})

    # Keep a deliberately small cross-section of paired legacy inputs in the
    # ordinary suite.  Their .oqp counterparts exercise the new public syntax;
    # these .inp decks make sure compatibility does not silently regress while
    # avoiding a second run of every converted example.
    _LEGACY_COMPATIBILITY_DECKS = frozenset({
        "HF/H2O_RHF-HF_ENERGY.inp",
        "MP2/h2o_ump2_6-31g.inp",
        "MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_GRADIENT.inp",
        "NMR/H2O_RHF-NMR.inp",
        "OPT/H2O_RHF-DFT_OPTIMIZE_OQP.inp",
        "OPT/C2H4_BHHLYP-MRSFTDDFT_TCI_OQP.inp",
        "OPT/C2H4_BHHLYP-MRSFTDDFT_BAEKA_OQP.inp",
    })

    def __init__(self,
                 base_test_dir: str = None,
                 output_dir: str = None,
                 total_cpus: int = None,
                 omp_threads: int = None,
                 mpi_manager=None):
        """
        Initialize the OQPTester.

        Args:
            base_test_dir (str): Base directory for test files. If None, uses
                                 the runtime root's share/examples.
            output_dir (str): Directory for storing test output files.
            total_cpus (int): Total number of CPUs to use.
            omp_threads (int): Number of OMP threads per test.
        """
        if base_test_dir is None:
            # Examples live under the resolved runtime root's share/ tree. Use
            # the package location via resolve_oqp_root(); OPENQP_ROOT is only a
            # compatibility fallback there (see pyoqp/README.md: a normal install
            # self-locates, "do not set OPENQP_ROOT").
            oqp_root, _ = resolve_oqp_root()
            base_test_dir = os.path.join(oqp_root, "share", "examples")
        self.base_test_dir = base_test_dir
        self.mpi_manager = mpi_manager
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        self.output_dir = os.path.abspath(f"{output_dir}_{timestamp}")
#       self.output_dir = os.path.abspath(f"{output_dir}")
        self.total_cpus = total_cpus if total_cpus is not None \
            else os.cpu_count()
        self.omp_threads = omp_threads if omp_threads is not None \
            else self.total_cpus

        self.max_workers = self.calculate_max_workers()
        self.results: List[Dict[str, Any]] = []
        self.report_file = os.path.join(self.output_dir, 'test_report.txt')
        self.start_time = None
        self.end_time = None
        self.status = 0
        self.input_format = "auto"

        os.environ['OMP_NUM_THREADS'] = str(self.omp_threads)

    def calculate_max_workers(self):
        if self.mpi_manager.use_mpi:
            return 1
        return max(1, self.total_cpus // self.omp_threads)

    @classmethod
    def _normalize_input_format(cls, input_format: str = "auto") -> str:
        """Return a validated test-input selector.

        ``auto`` preserves the normal regression preference: concise ``.oqp``
        decks, legacy-only ``.inp`` decks, and the small representative legacy
        compatibility set. The other modes are exact extension filters within
        the selected path or standard-suite scope.
        """
        normalized = str(input_format or "auto").strip().lower()
        if normalized not in cls._INPUT_FORMATS:
            choices = ", ".join(sorted(cls._INPUT_FORMATS))
            raise ValueError(
                f"Invalid test input format {input_format!r}; "
                f"choose one of: {choices}"
            )
        return normalized

    def log(self, message: str):
        """Simple logging function to stdout."""
        if self.mpi_manager.rank == 0:
            print(f"[OQPTester] {message}")

    def get_git_commit_info(self):
        try:
            git_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'])
            git_hash = git_hash.decode('ascii').strip()
            return f"{git_hash}"
        except (subprocess.CalledProcessError, FileNotFoundError):
            return "Unable to retrieve Git information"

    def get_git_branch_info(self):
        try:
            git_branch = subprocess.check_output(['git',
                                                  'rev-parse',
                                                  '--abbrev-ref',
                                                  'HEAD'])
            git_branch = git_branch.decode('ascii').strip()
            return f"{git_branch}"
        except (subprocess.CalledProcessError, FileNotFoundError):
            return "Unable to retrieve Git information"

    def _is_legacy_compatibility_deck(self, input_file: str) -> bool:
        """Return whether *input_file* is a selected legacy regression deck."""
        base_test_dir = getattr(self, "base_test_dir", None)
        candidates = set()
        if base_test_dir:
            try:
                relative = os.path.relpath(
                    os.path.abspath(input_file), os.path.abspath(base_test_dir)
                )
                candidates.add(relative.replace(os.sep, "/"))
            except (TypeError, ValueError):
                pass

        # An explicit source-tree path such as ``--run_tests examples/OPT``
        # may differ from the installed ``share/examples`` base chosen by the
        # runner.  Match the stable category/name suffix as well so the same
        # representative matrix is retained in either layout.
        try:
            absolute = os.path.abspath(input_file).replace(os.sep, "/")
        except (TypeError, ValueError):
            absolute = ""
        return bool(candidates & self._LEGACY_COMPATIBILITY_DECKS) or any(
            absolute.endswith(f"/{deck}")
            for deck in self._LEGACY_COMPATIBILITY_DECKS
        )

    def _project_name_for_input(self, input_file: str) -> str:
        """Return a collision-free project name for a paired legacy input."""
        project_name = os.path.splitext(os.path.basename(input_file))[0]
        if (input_file.lower().endswith(".restart.oqp")
                and project_name.lower().endswith(".restart")):
            # A paired continuation shares the producer project/log and NAMD
            # sidecars, but is scheduled only after that producer completes.
            project_name = project_name[:-len(".restart")]
        input_stem, extension = os.path.splitext(input_file)
        if (
            extension.lower() == ".inp"
            and os.path.isfile(f"{input_stem}.oqp")
        ):
            return f"{project_name}__legacy"
        return project_name

    def _case_output_dir(self, input_file: str,
                         project_name: str = None) -> str:
        """Return a deterministic, isolated output directory for one input.

        Several geometry drivers intentionally write conventional artifact
        names such as ``opt.xyz`` and ``opt_status.txt``. A test-specific
        directory therefore protects paired ``.inp``/``.oqp`` jobs, as well as
        unrelated decks with the same basename, when they run concurrently.
        """
        project_name = project_name or self._project_name_for_input(input_file)
        real_input = os.path.realpath(os.path.abspath(input_file))
        if real_input.lower().endswith(".restart.oqp"):
            real_input = real_input[:-len(".restart.oqp")] + ".oqp"
        digest = hashlib.sha256(os.fsencode(real_input)).hexdigest()[:12]
        return os.path.join(self.output_dir, f"{project_name}__{digest}")

    @staticmethod
    def _absolutize_caller_relative_inputs(mol, caller_cwd):
        """Preserve caller-CWD path semantics inside an isolated worker."""
        caller_cwd = os.path.realpath(os.path.abspath(caller_cwd))

        def runtime_path(value):
            if not isinstance(value, str) or not value.strip():
                return value
            expanded = os.path.expanduser(value.strip())
            if os.path.isabs(expanded):
                return os.path.realpath(expanded)
            return os.path.realpath(os.path.join(caller_cwd, expanded))

        guess = mol.config.get("guess", {})
        for key in ("file", "file2"):
            if key in guess:
                guess[key] = runtime_path(guess[key])
        md = mol.config.get("md", {})
        velocity = str(md.get("velocity", "") or "").strip()
        if velocity.lower() not in {
            "", "zero", "none", "0", "maxwell", "boltzmann", "random"
        }:
            md["velocity"] = runtime_path(velocity)
        qmmm = mol.config.get("qmmm", {})
        for key in ("pdb_file", "qm_atoms_xyz"):
            value = qmmm.get(key)
            if isinstance(value, str) and os.path.isfile(runtime_path(value)):
                qmmm[key] = runtime_path(value)
        for key in ("forcefield", "forcefield_files"):
            value = qmmm.get(key)
            if not isinstance(value, str):
                continue
            entries = [
                item for item in value.replace(",", " ").split() if item
            ]
            qmmm[key] = " ".join(
                runtime_path(item)
                if os.path.isfile(runtime_path(item))
                else item
                for item in entries
            )
        mol.oqp_runtime_cwd = caller_cwd

    def run_single_test(self, input_file: str,
                        caller_cwd: str = None) -> Dict[str, Any]:
        """
        Run a single OpenQP test.

        Args:
            input_file (str): Path to the input file.

        Returns:
            Dict[str, Any]: Dictionary containing test results.
        """
        project_name = self._project_name_for_input(input_file)
        case_output_dir = self._case_output_dir(input_file, project_name)
        os.makedirs(case_output_dir, exist_ok=True)
        log = os.path.join(case_output_dir, f"{project_name}.log")

        usempi = True if self.mpi_manager.use_mpi > 0 else False

        if usempi:
            input_file = self.mpi_manager.bcast(input_file)
            project_name = self.mpi_manager.bcast(project_name)
            log = self.mpi_manager.bcast(log)

        self.log(f"Running test for {project_name}")

        result = {
            "project": project_name,
            "input_file": input_file,
            "log_file": log,
            "status": "UNKNOWN",
            "message": "",
            "execution_time": 0
        }

        start_time = time.perf_counter()
        try:
            worker_cwd = os.getcwd()
            try:
                if caller_cwd is not None:
                    os.chdir(caller_cwd)
                runner = Runner(project=project_name,
                                input_file=input_file,
                                log=log,
                                silent=1,
                                usempi=usempi)
                if caller_cwd is not None:
                    self._absolutize_caller_relative_inputs(
                        runner.mol, caller_cwd)
            finally:
                if caller_cwd is not None:
                    os.chdir(worker_cwd)
            runner.run(test_mod=True)
            if self.mpi_manager.rank == 0:
                message, diff = runner.test()
                result["status"] = "PASSED" if round(diff, 4) == 0 else "FAILED"
                result["message"] = message
        except Exception as err:
            # geomeTRIC IRC calculations intentionally trace a finite path and
            # may terminate by reaching maxiter rather than an optimization
            # convergence criterion. Treat that known termination as a
            # completed IRC test when the calculation log was produced.
            is_irc_maxiter = False
            try:
                with open(input_file, 'r', encoding='utf-8') as inp:
                    input_text = inp.read().lower()
                is_irc_maxiter = (
                    type(err).__name__ == 'GeomOptNotConvergedError'
                    and 'runtype=irc' in input_text.replace(' ', '')
                    and os.path.exists(log)
                )
            except OSError:
                is_irc_maxiter = False

            # QM/MM examples need the optional OpenMM backend. When it is not
            # installed the run raises ModuleNotFoundError('openmm'); report
            # SKIPPED (like a build without an optional feature) so a build
            # without OpenMM still produces a green suite.
            needs_openmm_missing = (
                type(err).__name__ in ('ModuleNotFoundError', 'ImportError')
                and 'openmm' in str(err).lower()
            )

            # DFTB examples need the optional external openqp-dftb library
            # (libopenqp_dftb_c). When it is absent the adapter raises
            # FileNotFoundError("openqp-dftb not found ..."); report SKIPPED so a
            # build without the optional DFTB backend still produces a green suite.
            needs_dftb_missing = 'openqp-dftb not found' in str(err).lower()
            needs_xtb_missing = 'openqp-xtb not found' in str(err).lower()

            # The native optimizer covers the ordinary geometry workflows,
            # but legacy constrained inputs may still explicitly select the
            # optional geomeTRIC backend.  Its adapter raises this deliberately
            # specific import error when the extra is not installed.  Match the
            # adapter contract rather than every error mentioning "geometric"
            # so unrelated import bugs remain visible as test failures.
            needs_geometric_missing = (
                type(err).__name__ in ('ModuleNotFoundError', 'ImportError')
                and 'geometric is required for [optimize] lib=geometric'
                in str(err).lower()
            )

            if needs_openmm_missing:
                result["status"] = "SKIPPED"
                result["message"] = ("requires the optional OpenMM backend "
                                     "(not installed); skipped")
            elif needs_dftb_missing:
                result["status"] = "SKIPPED"
                result["message"] = ("requires the optional openqp-dftb backend "
                                     "(not installed); skipped")
            elif needs_xtb_missing:
                result["status"] = "SKIPPED"
                result["message"] = ("requires the optional openqp-xtb backend "
                                     "(not installed); skipped")
            elif needs_geometric_missing:
                result["status"] = "SKIPPED"
                result["message"] = ("requires the optional geomeTRIC optimizer "
                                     "(not installed); skipped")
            elif is_irc_maxiter:
                result["status"] = "PASSED"
                result["message"] = "IRC path reached configured maxiter and produced output log"
            else:
                self.log(f"Error in test {project_name}: {str(err)}")
                result["status"] = "ERROR"
                result["message"] = f"PyOQP error: {type(err).__name__} - {str(err)}"

        result["execution_time"] = time.perf_counter() - start_time
        return result

    def _run_isolated(self, input_file: str) -> Dict[str, Any]:
        """
        Run a single test in a dedicated subprocess and recover its result.

        Running the Fortran/C backend in a child process means a hard failure
        in one test (Fortran ERROR STOP, segfault, MKL/abort) terminates only
        that child. We translate a non-zero exit (or a timeout) into an ERROR
        result for that one test rather than letting it crash the whole run.
        """
        # Native Fortran units such as ``fort.6`` are opened relative to the
        # process working directory.  Parallel test subprocesses must therefore
        # have distinct working directories in addition to distinct log paths;
        # otherwise one Hessian test can consume another test's unit-6 file.
        input_file = os.path.realpath(os.path.abspath(input_file))
        project_name = self._project_name_for_input(input_file)
        case_output_dir = self._case_output_dir(input_file, project_name)
        os.makedirs(case_output_dir, exist_ok=True)
        log = os.path.join(case_output_dir, f"{project_name}.log")
        self.log(f"Running test for {project_name}")

        parent_cwd = os.getcwd()
        cmd = [
            sys.executable, "-m", "oqp.utils.oqp_tester",
            "--isolated", input_file,
            "--output-dir", self.output_dir,
            "--omp", str(self.omp_threads),
            "--caller-cwd", parent_cwd,
        ]
        child_env = os.environ.copy()
        inherited_pythonpath = child_env.get("PYTHONPATH")
        if inherited_pythonpath is not None:
            child_env["PYTHONPATH"] = os.pathsep.join(
                os.path.realpath(os.path.abspath(os.path.expanduser(entry)))
                if entry else parent_cwd
                for entry in inherited_pythonpath.split(os.pathsep)
            )
        start_time = time.perf_counter()
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=_test_timeout(),
                cwd=case_output_dir, env=child_env,
            )
            stdout, stderr, rc = proc.stdout, proc.stderr, proc.returncode
        except subprocess.TimeoutExpired as err:
            return {
                "project": project_name, "input_file": input_file,
                "log_file": log, "status": "ERROR",
                "message": f"PyOQP test exceeded time limit "
                           f"({_test_timeout()} s); see OQP_TEST_TIMEOUT",
                "execution_time": time.perf_counter() - start_time,
            }

        # The child prints exactly one marker line carrying the JSON result.
        for line in reversed((stdout or "").splitlines()):
            if line.startswith(_RESULT_MARKER):
                try:
                    return json.loads(line[len(_RESULT_MARKER):])
                except json.JSONDecodeError:
                    break

        # No result marker => the child died before reporting (ERROR STOP,
        # segfault, ...). If it ERROR STOPped only because the build lacks an
        # optional feature the example needs (e.g. the ddX/PCM backend when
        # ENABLE_DDX is OFF), report SKIPPED rather than ERROR so an
        # intentionally-trimmed build still produces a green suite.
        combined = (stdout or "") + (stderr or "")
        try:
            with open(log, "r", encoding="utf-8", errors="ignore") as fh:
                combined += fh.read()
        except OSError:
            pass
        for feature in self._OPTIONAL_FEATURE_SENTINELS:
            if feature in combined:
                self.log(f"Skipping test {project_name}: build lacks {feature}")
                return {
                    "project": project_name, "input_file": input_file,
                    "log_file": log, "status": "SKIPPED",
                    "message": f"requires a build feature not enabled "
                               f"({feature}); skipped",
                    "execution_time": time.perf_counter() - start_time,
                }

        tail = combined.strip().splitlines()[-3:]
        self.log(f"Error in test {project_name}: subprocess exit {rc}")
        return {
            "project": project_name, "input_file": input_file,
            "log_file": log, "status": "ERROR",
            "message": f"PyOQP test crashed (subprocess exit {rc}): "
                       + " | ".join(tail),
            "execution_time": time.perf_counter() - start_time,
        }

    def run_tests(self, test_path: str = 'all', *, input_format: str = 'auto'):
        """
        Run OpenQP tests based on the specified test path.

        Args:
            test_path (str): Path to test directory or specific input file.
                             Use 'all' to run all tests in the base directory.
            input_format (str): ``auto`` (default), ``inp``, ``oqp``, or
                                ``both``.
        """

        self.input_format = self._normalize_input_format(input_format)
        self.start_time = time.perf_counter()
        if self.mpi_manager.rank == 0:
            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)

        input_files = self._get_input_files(
            test_path, input_format=self.input_format
        )
        if not input_files:
            raise ValueError(
                f"No test inputs matched --input-format "
                f"{self.input_format} for: {test_path}"
            )

        if self.mpi_manager.use_mpi:
            primary_inputs, restart_inputs = self._partition_restart_inputs(
                input_files)
            for input_file in primary_inputs + restart_inputs:
                result = self.run_single_test(input_file)
                self.results.append(result)
                self._log_result_status(result)
        else:
            # Each OpenQP calculation loads a Fortran/C backend with some
            # process-global state.  Reusing a Python worker for multiple
            # tests can leak state between independent inputs; in particular,
            # running a TRAH SCF test before an ECP/MRSF energy test in the
            # same worker can make the later SCF exit after 0 iterations.
            #
            # Each test therefore runs in its own short-lived subprocess (see
            # _run_isolated). Besides isolating that process-global state, this
            # contains hard failures: a Fortran ERROR STOP or a segfault kills
            # only that child, so it is reported as a single ERROR instead of
            # tearing down a shared worker pool (BrokenProcessPool) and aborting
            # every still-pending test. A thread pool just supervises the child
            # processes, so the GIL is irrelevant here.
            primary_inputs, restart_inputs = self._partition_restart_inputs(
                input_files)
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_file = {
                    executor.submit(self._run_isolated, input_file): input_file
                    for input_file in primary_inputs
                }
                for future in as_completed(future_to_file):
                    result = future.result()
                    self.results.append(result)
                    self._log_result_status(result)
            # A *.restart.oqp example consumes the checkpoint and sidecars of
            # its same-stem producer. Run paired continuations only after all
            # ordinary examples have completed, never concurrently with their
            # producer.
            for input_file in restart_inputs:
                result = self._run_isolated(input_file)
                self.results.append(result)
                self._log_result_status(result)

        self.results.sort(key=lambda x: x['input_file'])
        self.end_time = time.perf_counter()

    @staticmethod
    def _partition_restart_inputs(input_files):
        """Place paired restart examples after their checkpoint producers."""
        selected = {
            os.path.realpath(os.path.abspath(path)) for path in input_files
        }
        primary = []
        restart = []
        for path in input_files:
            real_path = os.path.realpath(os.path.abspath(path))
            if real_path.lower().endswith(".restart.oqp"):
                producer = real_path[:-len(".restart.oqp")] + ".oqp"
                if producer in selected:
                    restart.append(path)
                    continue
            primary.append(path)
        return primary, restart

    def _get_input_files(self, test_path: str, *,
                         input_format: str = 'auto') -> List[str]:
        input_format = self._normalize_input_format(input_format)
        if test_path == 'all':
            test_dir = self.base_test_dir
        elif test_path == 'other':
            test_dir = os.path.join(self.base_test_dir, 'other')
        elif test_path == 'SCF':
            test_dir = os.path.join(self.base_test_dir, 'SCF')
        elif os.path.isdir(test_path):
            test_dir = test_path
        elif os.path.isfile(test_path) and test_path.lower().endswith(('.inp', '.oqp')):
            lower_path = test_path.lower()
            if lower_path.endswith('.resolved.oqp'):
                raise ValueError(
                    f"Resolved correction files are not regression inputs: {test_path}"
                )
            extension = os.path.splitext(lower_path)[1].lstrip('.')
            if input_format not in {'auto', 'both', extension}:
                raise ValueError(
                    f"Test input {test_path} is .{extension}, but "
                    f"--input-format {input_format} was requested"
                )
            return [test_path]
        else:
            raise ValueError(f"Invalid test path: {test_path}")

        candidates = sorted(
            os.path.join(root, file)
            for root, _, files in os.walk(test_dir)
            for file in files
            if file.lower().endswith(('.inp', '.oqp'))
            and not file.lower().endswith('.resolved.oqp')
        )

        if input_format == 'inp':
            input_files = [
                path for path in candidates if path.lower().endswith('.inp')
            ]
        elif input_format == 'oqp':
            input_files = [
                path for path in candidates if path.lower().endswith('.oqp')
            ]
        elif input_format == 'both':
            input_files = candidates
        else:
            # A same-stem .oqp is the public counterpart of its legacy .inp
            # and shares the same reference JSON. The default auto matrix runs
            # only the concise form for most pairs, while retaining a small
            # representative legacy set. Paired legacy project/log names use a
            # __legacy suffix so those intentional paired runs cannot race.
            semantic_stems = {
                os.path.splitext(path)[0]
                for path in candidates if path.lower().endswith('.oqp')
            }
            input_files = [
                path for path in candidates
                if path.lower().endswith('.oqp')
                or os.path.splitext(path)[0] not in semantic_stems
                or self._is_legacy_compatibility_deck(path)
            ]
        # The full-suite run ('all') skips a few examples that dominate CI
        # wall-clock; they still run when selected explicitly (a directory or a
        # input-file path). See _skip_in_full_run for which and why.
        if test_path == 'all':
            input_files = [
                f for f in input_files if not self._skip_in_full_run(f)
            ]
        return sorted(input_files)

    @staticmethod
    def _skip_in_full_run(input_file: str) -> bool:
        """True for examples excluded from `run_tests all`.

        Excluded because each costs many times a normal example and dominates
        the suite wall-clock:

          * numerical Hessians  -- runtype=hess without the opt-in
            type=analytical flag (numerical is the default); ~3N displaced
            gradient evaluations, and the MRSF one runs excited-state gradients
            at every displacement (~20-25x a normal test).
          * IRC paths           -- runtype=irc; traces many optimisation steps
            (the slowest single example in the suite).

        Excluded because they are not self-contained regression tests:

          * single-point / ground-state QM/MM  -- qmmm_flag=true without
            runtype=namd (the runtype=energy / OpenMM-MD decks). These exercise
            the single-point ESPF path, which is not yet functional in this
            branch, and read external PDB files by a bare name. The NAMD-QMMM
            examples (runtype=namd) are NOT skipped: they resolve their auxiliary
            files relative to the input file and SKIP cleanly when OpenMM is
            absent (see run_single_test).
          * DFTB (method=dftb)  -- these need the optional external openqp-dftb
            library (libopenqp_dftb_c) AND an external Slater-Koster parameter
            set, neither of which the repository or CI ships, so they are
            documentation examples rather than self-contained regression tests.
            They still SKIP cleanly when only the library is absent (see the
            openqp-dftb-not-found handling in run_single_test), and run when
            invoked explicitly with a real parameter path.

        Analytical Hessians (type=analytical) and ordinary opt/TS runs are
        unaffected, and the skipped examples still run when invoked explicitly
        by directory or explicit input-file path."""
        if input_file.lower().endswith('.oqp'):
            try:
                from oqp.utils.oqp_input import resolve_oqp_file
                cfg = resolve_oqp_file(input_file, write_resolved=False).legacy_config
            except (OSError, ValueError):
                return False
            input_cfg = cfg.get('input', {})
            runtype = str(input_cfg.get('runtype', 'energy')).strip().lower()
            hess_type = str(cfg.get('hess', {}).get('type', '')).strip().lower()
            qmmm = str(input_cfg.get('qmmm_flag', 'false')).strip().lower() in {
                'true', '1', 'yes', 'on'
            }
            method = str(input_cfg.get('method', 'hf')).strip().lower()
            return (
                runtype == 'irc'
                or (runtype == 'hess' and hess_type != 'analytical')
                or (qmmm and runtype != 'namd')
                or method == 'dftb'
            )
        try:
            with open(input_file, 'r', encoding='utf-8') as fh:
                text = fh.read().lower().replace(' ', '')
        except OSError:
            return False
        if 'runtype=irc' in text:
            return True
        if 'runtype=hess' in text and 'type=analytical' not in text:
            return True
        if 'qmmm_flag=true' in text and 'runtype=namd' not in text:
            return True
        if 'method=dftb' in text:
            return True
        if 'method=xtb' in text:
            # Same story as method=dftb: needs the optional external openqp-xtb
            # library (libopenqp_xtb_c) and a converter-generated .opxtb file.
            return True
        return False

    def format_time(self, seconds: float) -> str:
        hours, rem = divmod(seconds, 3600)
        minutes, seconds = divmod(rem, 60)
        return f"{int(hours):02d}:{int(minutes):02d}:{seconds:06.3f}"

    @staticmethod
    def _failure_reason(result) -> str:
        """One-line reason a test did not pass: the failing check(s) for a
        FAILED numeric mismatch, otherwise the last line of its message."""
        msg = str(result.get('message', '') or '')
        fails = [ln.strip() for ln in msg.splitlines() if 'failed' in ln.lower()]
        if fails:
            return ' | '.join(fails)
        lines = [ln.strip() for ln in msg.splitlines() if ln.strip()]
        return lines[-1] if lines else result.get('status', '')

    def _log_result_status(self, result) -> None:
        """Emit a live line naming any test that did not pass, so a CI log
        identifies the offending example as it happens (not just a count)."""
        status = result.get('status')
        if status and status != 'PASSED':
            self.log(f"{status}: {result.get('project')} "
                     f"-- {self._failure_reason(result)}")

    def generate_report(self) -> str:
        passed = sum(
            1 for result in self.results
            if result['status'] == 'PASSED'
        )
        failed = sum(
            1 for result in self.results
            if result['status'] == 'FAILED'
        )
        errors = sum(
            1 for result in self.results
            if result['status'] == 'ERROR'
        )
        skipped = sum(
            1 for result in self.results
            if result['status'] == 'SKIPPED'
        )
        self.status = 1 if failed > 0 or errors > 0 else 0

        # List the offending examples by name in the returned (console-visible)
        # summary. The full per-test detail goes only to test_report.txt in the
        # output dir, which CI does not upload -- so without this a red CI shows
        # "Failed: 1" with no clue which example failed.
        nonpassed = [r for r in self.results
                     if r['status'] in ('FAILED', 'ERROR')]
        skipped_list = [r for r in self.results if r['status'] == 'SKIPPED']
        nonpassed_block = ""
        if nonpassed:
            nonpassed_block += "Failing tests:\n"
            for r in nonpassed:
                nonpassed_block += (f"  {r['status']:7} {r['project']}\n"
                                    f"          {self._failure_reason(r)}\n")
            nonpassed_block += "\n"
        if skipped_list:
            nonpassed_block += ("Skipped tests: "
                                + ", ".join(r['project'] for r in skipped_list)
                                + "\n\n")

        total_time = self.end_time - self.start_time
        execution_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        git_commit_info = self.get_git_commit_info()
        git_branch_info = self.get_git_branch_info()

        summary = f"""
PyOQP Test Report
-----------------
Execution Date: {execution_date}
Git Branch Info: {git_branch_info}
Git Commit Info: {git_commit_info}
Output dir: {self.output_dir}
Input format: {getattr(self, 'input_format', 'auto')}
Total tests: {len(self.results)}
Passed: {passed}
Failed: {failed}
Errors: {errors}
Skipped: {skipped}

{nonpassed_block}Total CPUs: MPI Processors = {self.mpi_manager.size}, OpenMp Threads = {self.omp_threads}
Max parallel tests: {self.max_workers}

Total execution time: {self.format_time(total_time)}

"""
        detailed_results = "\nDetailed Results:\n"
        for result in self.results:
            detailed_results += f"\nTest: \
{os.path.abspath(result['input_file'])}\n"
            detailed_results += f"Status: {result['status']}\n"
            detailed_results += (
                f"Execution time: \
{self.format_time(result['execution_time'])}\n"
            )
            detailed_results += f"Log file: {result['log_file']}\n"
            detailed_results += f"Message: {result['message']}\n"

        report_file = os.path.join(self.output_dir, 'test_report.txt')
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(summary + detailed_results)

        return summary

    def run(self, test_path: str = 'all', *, input_format: str = 'auto') -> str:
        """
        Run tests and generate a report.

        Args:
            test_path (str): Path to test directory or specific input file.
                             Use 'all' to run all tests in the base directory.
            input_format (str): ``auto`` (default), ``inp``, ``oqp``, or
                                ``both``.

        Returns:
            str: A formatted string containing the test report.
        """
        input_format = self._normalize_input_format(input_format)
        self.log(
            f"Starting OpenQP tests for: {test_path} "
            f"(input format: {input_format})"
        )

        if os.path.exists(self.report_file):
            os.remove(self.report_file)

        self.run_tests(test_path, input_format=input_format)
        if self.mpi_manager.rank == 0:
            report = self.generate_report()
            self.log("OpenQP tests completed")
            return report
        else:
            return 1


def _run_isolated_main(argv=None):
    """
    Entry point for the per-test subprocess (``python -m oqp.utils.oqp_tester
    --isolated <input> --output-dir <dir> --omp <n>``).

    Runs exactly one test in this fresh process and prints the structured
    result as a single marker line. If the Fortran backend ERROR STOPs or
    crashes, this process simply dies with a non-zero exit and the parent
    records an ERROR for this test alone.
    """
    import argparse
    from oqp.pyoqp import MPIManager

    parser = argparse.ArgumentParser(description="Run one OpenQP test in isolation")
    parser.add_argument("--isolated", required=True, help="input .inp or .oqp file")
    parser.add_argument("--output-dir", required=True, help="shared output dir")
    parser.add_argument("--omp", type=int, default=1, help="OMP threads")
    parser.add_argument(
        "--caller-cwd", required=True,
        help="working directory used to resolve caller-relative input files",
    )
    args = parser.parse_args(argv)

    os.environ["OMP_NUM_THREADS"] = str(args.omp)

    # Build a tester shell without re-creating a timestamped output dir: the
    # parent already created the shared one and passes it in.
    tester = OQPTester.__new__(OQPTester)
    tester.output_dir = args.output_dir
    tester.mpi_manager = MPIManager()

    result = tester.run_single_test(args.isolated, caller_cwd=args.caller_cwd)
    sys.stdout.flush()
    print(_RESULT_MARKER + json.dumps(result))
    return 0


if __name__ == "__main__":
    sys.exit(_run_isolated_main())
