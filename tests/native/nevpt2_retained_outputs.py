"""Oracle for original native outputs retained across later calls and errors."""
import numpy as np


class RetainedOutputs:
    def __init__(self):
        self._records = []

    def __len__(self):
        return len(self._records)

    def remember(self, outputs):
        # Keep the actual arrays (including borrowed views), alongside owned snapshots.
        original = tuple(outputs)
        assert all(np.all(np.isfinite(output)) for output in original)
        self._records.append((original, tuple(np.array(x, copy=True) for x in original)))

    def assert_unchanged(self):
        for original, snapshots in self._records:
            for actual, expected in zip(original, snapshots):
                np.testing.assert_array_equal(actual, expected)
