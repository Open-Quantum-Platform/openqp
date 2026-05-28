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
import time
import subprocess
import shutil
from typing import List, Dict, Any
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

from oqp.pyoqp import Runner

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
    def __init__(self,
                 base_test_dir: str = None,
                 output_dir: str = None,
                 total_cpus: int = None,
                 omp_threads: int = None,
                 mpi_manager=None):
        """
        Initialize the OQPTester.

        Args:
            base_test_dir (str): Base directory for test files.
                                 If None, uses $OPENQP_ROOT/examples.
            output_dir (str): Directory for storing test output files.
            total_cpus (int): Total number of CPUs to use.
            omp_threads (int): Number of OMP threads per test.
        """
        try:
            oqp_root = os.environ["OPENQP_ROOT"]
        except KeyError:
            oqp_root =  os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
        if not oqp_root:
            raise EnvironmentError("OPENQP_ROOT environment variable is not set")

        self.base_test_dir = base_test_dir or os.path.join(oqp_root,
                                                           'share/examples')
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

    @staticmethod
    def _read_simple_input_value(input_file: str, section: str, key: str, default: str = "") -> str:
        current_section = None
        try:
            with open(input_file, encoding="utf-8") as handle:
                for raw_line in handle:
                    line = raw_line.strip()
                    if not line or line.startswith(("#", ";")):
                        continue
                    if line.startswith("[") and line.endswith("]"):
                        current_section = line[1:-1].strip().lower()
                        continue
                    if current_section == section.lower() and "=" in line:
                        name, value = line.split("=", 1)
                        if name.strip().lower() == key.lower():
                            return value.strip().strip('"').strip("'")
        except OSError:
            return default
        return default

    def _optional_backend_skip_reason(self, input_file: str) -> str:
        method = self._read_simple_input_value(input_file, "input", "method").lower()
        if method != "dftb":
            return ""
        executable = self._read_simple_input_value(input_file, "dftb", "executable", "dftb+") or "dftb+"
        if shutil.which(executable) is None:
            return (
                f"DFTB+ executable not found on PATH: {executable}. "
                "Skipping optional DFTB+ backend test; install DFTB+ and configure "
                "the [dftb] executable/sk_path settings for full test coverage."
            )

        sk_path = self._read_simple_input_value(input_file, "dftb", "sk_path")
        sk_path = os.path.expandvars(os.path.expanduser(sk_path))
        if not sk_path:
            return (
                "DFTB+ parameter directory is not configured; set [dftb] sk_path. "
                "Skipping optional DFTB+ backend test; configure the Slater-Koster "
                "parameter directory for full test coverage."
            )
        if not os.path.isdir(sk_path):
            return (
                f"DFTB+ parameter directory not found: {sk_path}. "
                "Skipping optional DFTB+ backend test; configure the [dftb] "
                "sk_path setting for full test coverage."
            )
        return ""

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
        skip_reason = self._optional_backend_skip_reason(input_file)
        if skip_reason:
            self.log(f"Skipping test {project_name}: {skip_reason}")
            result["status"] = "SKIPPED"
            result["message"] = skip_reason
            result["execution_time"] = time.perf_counter() - start_time
            return result

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
            self.log(f"Error in test {project_name}: {str(err)}")
            result["status"] = "ERROR"
            result["message"] = f"PyOQP error: {type(err).__name__} - {str(err)}"

        result["execution_time"] = time.perf_counter() - start_time
        return result

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
        if test_path == 'all':
            self.max_workers = int(os.environ.get('OQP_TEST_MAX_WORKERS', '1'))
        if not input_files:
            self.end_time = time.perf_counter()
            return

        if self.mpi_manager.use_mpi:
            for input_file in input_files:
                result = self.run_single_test(input_file)
                self.results.append(result)
        else:
            with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_file = {
                    executor.submit(self.run_single_test, input_file): input_file
                    for input_file in input_files
                }
                for future in as_completed(future_to_file):
                    input_file = future_to_file[future]
                    try:
                        self.results.append(future.result())
                    except Exception as err:
                        project_name = os.path.splitext(os.path.basename(input_file))[0]
                        self.log(f"Error in test {project_name}: {type(err).__name__} - {str(err)}")
                        self.results.append({
                            "project": project_name,
                            "input_file": input_file,
                            "log_file": os.path.join(self.output_dir, f"{project_name}.log"),
                            "status": "ERROR",
                            "message": f"PyOQP worker error: {type(err).__name__} - {str(err)}",
                            "execution_time": 0,
                        })

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

        return [
            os.path.join(root, file)
            for root, _, files in os.walk(test_dir)
            for file in files
            if file.endswith('.inp')
        ]

    def format_time(self, seconds: float) -> str:
        hours, rem = divmod(seconds, 3600)
        minutes, seconds = divmod(rem, 60)
        return f"{int(hours):02d}:{int(minutes):02d}:{seconds:06.3f}"

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

Total CPUs: MPI Processors = {self.mpi_manager.size}, OpenMp Threads = {self.omp_threads}
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

        try:
            self.run_tests(test_path)
        except Exception as err:
            self.log(f"OpenQP test harness error: {type(err).__name__} - {str(err)}")
            self.results.append({
                "project": "test_harness",
                "input_file": test_path,
                "log_file": self.report_file,
                "status": "ERROR",
                "message": f"PyOQP test harness error: {type(err).__name__} - {str(err)}",
                "execution_time": 0,
            })
            if self.start_time is None:
                self.start_time = time.perf_counter()
            self.end_time = time.perf_counter()
        if self.mpi_manager.rank == 0:
            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)
            report = self.generate_report()
            self.log("OpenQP tests completed")
            return report
        else:
            return 1
