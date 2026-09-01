"""Regression guards for conservative MRSF response-integral screening."""

from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "source" / "tdhf_mrsf_lib.F90"


def _compact_source() -> str:
    return "".join(LIB.read_text().lower().split())


def _shell_envelope(density: np.ndarray, shells: tuple[slice, ...]) -> np.ndarray:
    """Reference shell maximum over vector, component, and AO indices."""
    envelope = np.zeros((len(shells), len(shells)))
    for i, aos_i in enumerate(shells):
        for j, aos_j in enumerate(shells[: i + 1]):
            value = np.max(np.abs(density[:, :, aos_j, aos_i]))
            envelope[i, j] = value
            envelope[j, i] = value
    return envelope


def test_screening_uses_the_complete_response_density_batch():
    source = _compact_source()
    routine = source.split("subroutineint2_mrsf_data_t_init_screen", 1)[1]
    routine = routine.split("endsubroutine", 1)[0]
    shell_routine = source.split("subroutineshell_den_screen_mrsf", 1)[1]
    shell_routine = shell_routine.split("endsubroutine", 1)[0]

    assert "callshell_den_screen_mrsf(this%dsh,this%d3,basis)" in routine
    assert "this%d3(:,sized,:,:)" not in routine
    assert "dimension(:,:,:,:)::da" in shell_routine
    assert "maxval(abs(da(:,:,minj:maxj,mini:maxi)))" in shell_routine


def test_nonfinal_component_and_single_vector_control_the_shell_envelope():
    density = np.zeros((2, 7, 3, 3))
    shells = (slice(0, 1), slice(1, 3))

    # The largest density is deliberately in the first vector and a
    # nonfinal spin-adapted component.  Looking only at component 7 would
    # discard the shell pair required by this response contraction.
    density[0, 2, 0, 2] = -4.25
    density[1, 6, 0, 2] = 0.5
    envelope = _shell_envelope(density, shells)

    assert envelope[1, 0] == 4.25
    assert envelope[0, 1] == 4.25


def test_shell_envelope_is_invariant_to_batch_partitioning():
    density = np.zeros((3, 7, 3, 3))
    shells = (slice(0, 1), slice(1, 3))
    density[0, 1, 0, 1] = 1.25
    density[1, 4, 0, 2] = -3.5
    density[2, 6, 2, 2] = 2.0

    full = _shell_envelope(density, shells)
    partitioned = np.maximum.reduce(
        [_shell_envelope(density[i : i + 1], shells) for i in range(3)]
    )

    np.testing.assert_array_equal(full, partitioned)
