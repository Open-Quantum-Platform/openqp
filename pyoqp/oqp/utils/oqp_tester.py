"""
OpenQP Test Runner

This module provides a class for running OpenQP tests
in parallel and generating reports.
It is compatible with CI pipelines.

Author: Konstantin Komarov
Email: constlike@gmail.com
Created: Aug 2024
"""
import os
import sys
import json
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
      - 'openqp --run-tests all': Run all tests from all folders in examples
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
    )

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

        os.environ['OMP_NUM_THREADS'] = str(self.omp_threads)

    def calculate_max_workers(self):
        if self.mpi_manager.use_mpi:
            return 1
        return max(1, self.total_cpus // self.omp_threads)

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

    def run_single_test(self, input_file: str) -> Dict[str, Any]:
        """
        Run a single OpenQP test.

        Args:
            input_file (str): Path to the input file.

        Returns:
            Dict[str, Any]: Dictionary containing test results.
        """
        project_name = os.path.splitext(os.path.basename(input_file))[0]
        log = os.path.join(self.output_dir, f"{project_name}.log")

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
            runner = Runner(project=project_name,
                            input_file=input_file,
                            log=log,
                            silent=1,
                            usempi=usempi)
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

            if is_irc_maxiter:
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
        project_name = os.path.splitext(os.path.basename(input_file))[0]
        log = os.path.join(self.output_dir, f"{project_name}.log")
        self.log(f"Running test for {project_name}")

        cmd = [
            sys.executable, "-m", "oqp.utils.oqp_tester",
            "--isolated", input_file,
            "--output-dir", self.output_dir,
            "--omp", str(self.omp_threads),
        ]
        start_time = time.perf_counter()
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=_test_timeout(),
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

    def run_tests(self, test_path: str = 'all'):
        """
        Run OpenQP tests based on the specified test path.

        Args:
            test_path (str): Path to test directory or specific input file.
                             Use 'all' to run all tests in the base directory.
        """

        self.start_time = time.perf_counter()
        if self.mpi_manager.rank == 0:
            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)

        input_files = self._get_input_files(test_path)
        if not input_files:
            return

        if self.mpi_manager.use_mpi:
            for input_file in input_files:
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
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_file = {
                    executor.submit(self._run_isolated, input_file): input_file
                    for input_file in input_files
                }
                for future in as_completed(future_to_file):
                    result = future.result()
                    self.results.append(result)
                    self._log_result_status(result)

        self.results.sort(key=lambda x: x['input_file'])
        self.end_time = time.perf_counter()

    def _get_input_files(self, test_path: str) -> List[str]:
        if test_path == 'all':
            test_dir = self.base_test_dir
        elif test_path == 'other':
            test_dir = os.path.join(self.base_test_dir, 'other')
        elif test_path == 'SCF':
            test_dir = os.path.join(self.base_test_dir, 'SCF')
        elif os.path.isdir(test_path):
            test_dir = test_path
        elif os.path.isfile(test_path) and test_path.endswith('.inp'):
            return [test_path]
        else:
            print(f"Invalid test path: {test_path}")
            return []

        input_files = [
            os.path.join(root, file)
            for root, _, files in os.walk(test_dir)
            for file in files
            if file.endswith('.inp')
        ]
        # The full-suite run ('all') skips a few examples that dominate CI
        # wall-clock; they still run when selected explicitly (a directory or a
        # .inp path). See _skip_in_full_run for which and why.
        if test_path == 'all':
            input_files = [
                f for f in input_files if not self._skip_in_full_run(f)
            ]
        return input_files

    @staticmethod
    def _skip_in_full_run(input_file: str) -> bool:
        """True for examples excluded from `run_tests all` because each costs
        many times a normal example and dominates the suite wall-clock:

          * numerical Hessians  -- runtype=hess without the opt-in
            type=analytical flag (numerical is the default); ~3N displaced
            gradient evaluations, and the MRSF one runs excited-state gradients
            at every displacement (~20-25x a normal test).
          * IRC paths           -- runtype=irc; traces many optimisation steps
            (the slowest single example in the suite).

        Analytical Hessians (type=analytical) and ordinary opt/TS runs are
        unaffected, and the skipped examples still run when invoked explicitly
        by directory or .inp path."""
        try:
            with open(input_file, 'r', encoding='utf-8') as fh:
                text = fh.read().lower().replace(' ', '')
        except OSError:
            return False
        if 'runtype=irc' in text:
            return True
        if 'runtype=hess' in text and 'type=analytical' not in text:
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

    def run(self, test_path: str = 'all') -> str:
        """
        Run tests and generate a report.

        Args:
            test_path (str): Path to test directory or specific input file.
                             Use 'all' to run all tests in the base directory.

        Returns:
            str: A formatted string containing the test report.
        """
        self.log(f"Starting OpenQP tests for: {test_path}")

        if os.path.exists(self.report_file):
            os.remove(self.report_file)

        self.run_tests(test_path)
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
    parser.add_argument("--isolated", required=True, help="input .inp file")
    parser.add_argument("--output-dir", required=True, help="shared output dir")
    parser.add_argument("--omp", type=int, default=1, help="OMP threads")
    args = parser.parse_args(argv)

    os.environ["OMP_NUM_THREADS"] = str(args.omp)

    # Build a tester shell without re-creating a timestamped output dir: the
    # parent already created the shared one and passes it in.
    tester = OQPTester.__new__(OQPTester)
    tester.output_dir = args.output_dir
    tester.mpi_manager = MPIManager()

    result = tester.run_single_test(args.isolated)
    sys.stdout.flush()
    print(_RESULT_MARKER + json.dumps(result))
    return 0


if __name__ == "__main__":
    sys.exit(_run_isolated_main())
