"""Exercise the actual Fock kernels and reduction, with an independent signed oracle.

Run with the clean installed package's Python, at OMP_NUM_THREADS=1 and 4.
For a compiled MPI package also run with mpiexec -n 2 and
OQP_SHARED_REQUIRE_MPI=1; a serial library fails that mode, rather than skipping.
"""
import os
import re
import unittest
from pathlib import Path


class SharedMRSFFockTests(unittest.TestCase):
    def test_single_thread_exchange_has_no_atomic_updates(self):
        """Reject serial synchronization without timing-dependent assertions.

        The signed native oracle separately executes both thread-count branches.
        Here every atomic loop must have an identical ordinary serial branch.
        """
        source = (Path(__file__).resolve().parents[1] /
                  "source" / "tdhf_mrsf_lib.F90").read_text()
        for name in ("int2_mrsf_update_worker", "int2_umrsf_update_worker",
                     "int2_mrsf_shell_worker", "int2_umrsf_shell_worker"):
            with self.subTest(worker=name):
                worker = re.search(r"(?ms)^  subroutine " + name +
                                   r"\(.*?^  end subroutine[^\n]*", source).group()
                dispatches = list(re.finditer(
                    r"(?ms)^([ \t]*)if \(this%nthreads == 1\) then\n"
                    r"(.*?)^\1else\n(.*?)^\1end if", worker))
                self.assertEqual(len(dispatches), 3,
                                 "all three exchange loops need serial dispatch")
                covered = 0
                for dispatch in dispatches:
                    ordinary, atomic = dispatch.group(2), dispatch.group(3)
                    self.assertNotIn("!$omp atomic", ordinary)
                    self.assertEqual(atomic.count("!$omp atomic update"), 8,
                                     "threaded exchange must remain synchronized")
                    covered += atomic.count("!$omp atomic update")
                    stripped = "\n".join(line for line in atomic.splitlines()
                                         if "!$omp atomic update" not in line)
                    self.assertEqual(ordinary.rstrip(), stripped.rstrip())
                self.assertEqual(covered, worker.count("!$omp atomic update"),
                                 "atomic exchange can run without serial dispatch")

    def test_private_storage_is_the_default(self):
        source = (Path(__file__).resolve().parents[1] /
                  "source" / "tdhf_mrsf_lib.F90").read_text()
        self.assertRegex(source, r"logical :: shared_exchange = \.false\.")
        for entry, private in (("int2_mrsf_data_t_update", "int2_mrsf_private_update"),
                               ("int2_umrsf_data_t_update", "int2_umrsf_private_update"),
                               ("int2_mrsf_consume_shell", "int2_mrsf_private_shell"),
                               ("int2_umrsf_consume_shell", "int2_umrsf_private_shell")):
            routine = re.search(r"(?ms)^  subroutine " + entry +
                                r"\(.*?^  end subroutine[^\n]*", source).group()
            self.assertIn("if (this%shared_exchange) then", routine)
            self.assertIn("call " + private + "(", routine)

    def test_compact_api_and_installed_control_select_storage(self):
        from oqp import OpenQP
        from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA, OQPData

        self.assertEqual(OQP_CONFIG_SCHEMA["tdhf"]["shared_f3"]["default"], "False")
        job = OpenQP(project="shared_f3_api").molecule("H 0 0 0; H 0 0 0.74")
        job.settings.tdhf(shared_f3=True)
        self.assertEqual(job.to_input_dict()["tdhf"]["shared_f3"], "True")
        data = OQPData(silent=1)
        self.assertFalse(data._data.control.mrsf_shared_f3)
        data.parse_section({"tdhf": {"shared_f3": True}}, "tdhf")
        self.assertTrue(data._data.control.mrsf_shared_f3)
        data.parse_section({"tdhf": {"shared_f3": False}}, "tdhf")
        self.assertFalse(data._data.control.mrsf_shared_f3)

    def test_signed_channels_fp32_cam_resize_and_cleanup(self):
        import oqp

        # Exercise the generated binding: this ABI takes no molecule handle.
        # A missing no-handle registration must fail here, not be bypassed.
        check = oqp.int2_mrsf_shared_selftest
        use_mpi = os.getenv("OQP_SHARED_REQUIRE_MPI") == "1"
        comm, ranks = 0, 1
        if use_mpi:
            from mpi4py import MPI
            comm, ranks = MPI.COMM_WORLD.py2f(), MPI.COMM_WORLD.size
            self.assertGreaterEqual(ranks, 2, "multirank gate needs at least two ranks")
        threads = int(os.environ["OMP_NUM_THREADS"])
        self.assertIn(threads, (1, 4))
        for _ in range(2):
            error = oqp.ffi.new("double *")
            checks, failures = oqp.ffi.new("int *"), oqp.ffi.new("int *")
            check(comm, use_mpi, ranks, threads, error, checks, failures)
            self.assertEqual(checks[0], 78)
            self.assertEqual(failures[0], 0)
            self.assertEqual(error[0], 0.0)


if __name__ == "__main__":
    unittest.main()
