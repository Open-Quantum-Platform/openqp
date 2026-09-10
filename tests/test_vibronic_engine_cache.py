"""HarmonicOverlapEngine memoization must not outlive the engine."""

import gc
import importlib.util
from pathlib import Path
import sys
import unittest
import weakref

ROOT = Path(__file__).resolve().parents[1]


def load_vibronic_module():
    path = ROOT / "pyoqp" / "oqp" / "library" / "vibronic.py"
    spec = importlib.util.spec_from_file_location("openqp_vibronic_engine_cache", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class EngineCacheTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vibronic = load_vibronic_module()

    def model(self):
        return self.vibronic.HarmonicVibronicModel.create(
            [1000.0],
            [900.0],
            [[1.0]],
            [0.3],
            coordinate_unit="sqrt(me)*bohr",
            coordinate_phase_convention="synthetic cache test",
        )

    def overlaps(self, engine):
        return [engine.overlap((i,), (j,)) for i in range(3) for j in range(4)]

    def test_spent_engines_are_released(self):
        # A class-level lru_cache on the method kept every engine (and every
        # polynomial it evaluated) alive through the cache key.
        engine = self.vibronic.HarmonicOverlapEngine(self.model())
        self.overlaps(engine)
        released = weakref.ref(engine)
        del engine
        gc.collect()
        self.assertIsNone(released(), "a spent HarmonicOverlapEngine is still referenced")

    def test_memoized_overlaps_are_unchanged(self):
        engine = self.vibronic.HarmonicOverlapEngine(self.model())
        first = self.overlaps(engine)
        self.assertEqual(first, self.overlaps(engine))
        self.assertEqual(first, self.overlaps(self.vibronic.HarmonicOverlapEngine(self.model())))


if __name__ == "__main__":
    unittest.main()
