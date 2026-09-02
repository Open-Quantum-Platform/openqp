"""Determinant-free reference for the two-SOMO MRSF seven-density action.

The reference contains only the spin-coupled coefficient table of SI Secs. 5
and 8 and spatial-orbital tensor contractions.  It neither imports OpenQP nor
constructs a many-electron basis.  A tiny Fortran executable compiles the exact
``mrsfcbc`` and ``mrsfmntoia`` routines extracted from the checkout at test
time, so the numerical comparison exercises production source without keeping
a second copy of it in the test suite.

The minimal model has one closed orbital C, two open orbitals O1/O2, one
virtual orbital V, and four electrons.  Its expanded packed dimension is nine;
the physical singlet and triplet dimensions are eight and six, respectively.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import platform
import re
import shutil
import subprocess

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
FORTRAN_DRIVER = (
    ROOT / "tests/fortran/test_mrsf_spin_adapted_sigma_oracle.F90"
)
MRSF_SOURCE = ROOT / "source/tdhf_mrsf_lib.F90"
CONVENTION_SOURCE = ROOT / "source/modules/tdhf_mrsf_conventions.F90"

C, O1, O2, V = range(4)
NBF = 4
NPACKED = 9
ISQRT2 = 1.0 / np.sqrt(2.0)

PACKED_SLOTS = {
    "CO1": 0,
    "LR": 1,
    "G": 2,
    "CO2": 3,
    "D": 4,
    "R": 5,
    "CV": 6,
    "OV1": 7,
    "OV2": 8,
}
PHYSICAL_LABELS = {
    1: ("CO1", "CO2", "CV", "OV1", "OV2", "G", "D", "LR"),
    3: ("CO1", "CO2", "CV", "OV1", "OV2", "LR"),
}
FAMILY_CASES = {"a0": 0, "coco": 1, "ovov": 2, "coov": 3}


def _dyad(row: int, column: int) -> np.ndarray:
    matrix = np.zeros((NBF, NBF))
    matrix[row, column] = 1.0
    return matrix


E = {
    (row, column): _dyad(row, column)
    for row in range(NBF)
    for column in range(NBF)
}


def _density_coefficients(
    multiplicity: int,
    *,
    mutations: frozenset[str] = frozenset(),
) -> np.ndarray:
    """Return the explicit seven-density coefficient table.

    Axis order is physical coordinate, density channel, AO row, AO column.
    The table is a direct spin-coupled CSF/diagrammatic statement, not a
    translation of the production DGEMM sequence.
    """

    labels = PHYSICAL_LABELS[multiplicity]
    result = np.zeros((len(labels), 7, NBF, NBF))
    for column, label in enumerate(labels):
        if label == "CO1":
            result[column, 2] = E[C, O1]
            result[column, 5] = E[O2, C]
        elif label == "CO2":
            result[column, 3] = E[C, O2]
            sign = 1.0 if "d6_plus" in mutations else -1.0
            result[column, 5] = sign * E[O1, C]
        elif label == "CV":
            result[column, 6] = E[C, V]
        elif label == "OV1":
            result[column, 1] = E[O1, V]
            sign = 1.0 if "d5_plus" in mutations else -1.0
            result[column, 4] = sign * E[V, O2]
        elif label == "OV2":
            result[column, 0] = E[O2, V]
            result[column, 4] = E[V, O1]
        elif label == "G":
            result[column, 6] = E[O2, O1]
        elif label == "D":
            result[column, 6] = E[O1, O2]
        elif label == "LR":
            scale = 1.0 if "lr_input_unnormalized" in mutations else ISQRT2
            if multiplicity == 1:
                sign = 1.0 if "lr_input_wrong_sign" in mutations else -1.0
            else:
                sign = -1.0 if "lr_input_wrong_sign" in mutations else 1.0
            result[column, 6] = scale * (E[O1, O1] + sign * E[O2, O2])
        else:  # pragma: no cover - closed explicit table
            raise AssertionError(label)

        if label in {"CO1", "CO2", "OV1", "OV2"}:
            result[column, 6] += result[column, 0]
            result[column, 6] += result[column, 1]
            result[column, 6] += result[column, 2]
            result[column, 6] += result[column, 3]

    if "swap_channels_1_2" in mutations:
        result[:, [0, 1]] = result[:, [1, 0]]
    if "swap_channels_3_4" in mutations:
        result[:, [2, 3]] = result[:, [3, 2]]
    return result


def _pair_index(first: int, second: int) -> int:
    target = tuple(sorted((first, second)))
    pairs = [(p, q) for p in range(NBF) for q in range(p, NBF)]
    return pairs.index(target) + 1


def _synthetic_eri() -> np.ndarray:
    """Return a fixed eightfold-symmetric spatial ERI tensor."""

    eri = np.empty((NBF, NBF, NBF, NBF))
    for p in range(NBF):
        for q in range(NBF):
            for r in range(NBF):
                for s in range(NBF):
                    left = _pair_index(p, q)
                    right = _pair_index(r, s)
                    low, high = sorted((left, right))
                    eri[p, q, r, s] = (
                        17 * low + 31 * high + 7 * low * high
                    ) / 997.0
    return eri


ERI = _synthetic_eri()


def _diagrammatic_images(densities: np.ndarray) -> np.ndarray:
    """Contract all seven densities without using production code."""

    transposed = densities.swapaxes(-1, -2)
    direct = np.einsum(
        "mnkl,vckl->vcmn", ERI, densities + transposed, optimize=True
    )
    exchange = np.einsum(
        "mknl,vckl->vcmn", ERI, densities, optimize=True
    )
    images = -exchange
    images[:, :4] += direct[:, :4]
    return images


def _select_family(
    multiplicity: int,
    images: np.ndarray,
    family: str,
) -> np.ndarray:
    selected = np.zeros_like(images)
    selected[:, 6] = images[:, 6]
    sign = 1.0 if multiplicity == 1 else -1.0
    if family == "coco":
        selected[:, 5] = sign * images[:, 5]
    elif family == "ovov":
        selected[:, 4] = sign * images[:, 4]
    elif family == "coov":
        selected[:, :4] = sign * images[:, :4]
    elif family != "a0":  # pragma: no cover - closed public cases
        raise ValueError(f"unknown channel family: {family}")
    return selected


def _backproject(
    multiplicity: int,
    images: np.ndarray,
    *,
    mutations: frozenset[str] = frozenset(),
) -> np.ndarray:
    """Project AO images onto the physical packed coordinates."""

    if "flip_channel_7" in mutations:
        images = images.copy()
        images[:, 6] *= -1.0

    nvector = images.shape[0]
    result = np.zeros((NPACKED, nvector))
    f1, f2, f3, f4, f5, f6, f7 = np.moveaxis(images, 1, 0)

    d6_sign = 1.0 if "d6_output_plus" in mutations else -1.0
    d5_sign = 1.0 if "d5_output_plus" in mutations else -1.0
    result[PACKED_SLOTS["CO1"]] = f7[:, C, O1] + f1[:, C, O1] + d6_sign * f6[:, C, O2]
    result[PACKED_SLOTS["CO2"]] = f7[:, C, O2] + f2[:, C, O2] + f6[:, C, O1]
    result[PACKED_SLOTS["CV"]] = f7[:, C, V]
    result[PACKED_SLOTS["OV1"]] = f7[:, O1, V] + f4[:, O1, V] + f5[:, O2, V]
    result[PACKED_SLOTS["OV2"]] = f7[:, O2, V] + f3[:, O2, V] + d5_sign * f5[:, O1, V]

    lr_scale = 1.0 if "lr_output_unnormalized" in mutations else ISQRT2
    if multiplicity == 1:
        result[PACKED_SLOTS["G"]] = f7[:, O2, O1]
        result[PACKED_SLOTS["D"]] = f7[:, O1, O2]
        sign = 1.0 if "lr_output_wrong_sign" in mutations else -1.0
    else:
        if "retain_triplet_gd" in mutations:
            result[PACKED_SLOTS["G"]] = f7[:, O2, O1]
            result[PACKED_SLOTS["D"]] = f7[:, O1, O2]
        sign = -1.0 if "lr_output_wrong_sign" in mutations else 1.0
    result[PACKED_SLOTS["LR"]] = lr_scale * (
        f7[:, O1, O1] + sign * f7[:, O2, O2]
    )
    return result


def _oracle_sigma(
    multiplicity: int,
    family: str,
    *,
    density_mutations: frozenset[str] = frozenset(),
    output_mutations: frozenset[str] = frozenset(),
) -> tuple[np.ndarray, np.ndarray]:
    densities = _density_coefficients(
        multiplicity, mutations=density_mutations
    )
    images = _diagrammatic_images(densities)
    selected = _select_family(multiplicity, images, family)
    sigma = _backproject(
        multiplicity, selected, mutations=output_mutations
    )
    return densities, sigma


def _physical_rows(multiplicity: int) -> tuple[int, ...]:
    return tuple(PACKED_SLOTS[label] for label in PHYSICAL_LABELS[multiplicity])


def _extract_subroutine(name: str) -> str:
    source = MRSF_SOURCE.read_text()
    match = re.search(
        rf"(?ms)^\s*subroutine\s+{name}\b.*?^\s*end\s+subroutine\s+{name}\s*$",
        source,
    )
    if match is None:
        raise AssertionError(f"production subroutine {name} was not found")
    return match.group(0)


def _production_extract_source() -> str:
    preamble = r"""
module precision
  implicit none
  integer, parameter :: dp=kind(1.0d0)
end module precision

module io_constants
  implicit none
  integer, parameter :: iw=6
end module io_constants

module messages
  implicit none
  integer, parameter :: with_abort=1
contains
  subroutine show_message(message,mode)
    character(len=*), intent(in) :: message
    integer, intent(in) :: mode
    write(*,*) trim(message),mode
    error stop
  end subroutine show_message
end module messages

module types
  implicit none
  type basis_information
    integer :: nbf=0
  end type basis_information
  type molecular_properties
    integer :: nelec_a=0,nelec_b=0
  end type molecular_properties
  type tddft_information
    integer :: mult=0
    logical :: debug_mode=.false.
  end type tddft_information
  type information
    type(basis_information) :: basis
    type(molecular_properties) :: mol_prop
    type(tddft_information) :: tddft
  end type information
end module types

subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  use precision, only: dp
  implicit none
  character(len=1), intent(in) :: transa,transb
  integer, intent(in) :: m,n,k,lda,ldb,ldc
  real(kind=dp), intent(in) :: alpha,beta
  real(kind=dp), intent(in) :: a(lda,*),b(ldb,*)
  real(kind=dp), intent(inout) :: c(ldc,*)
  real(kind=dp) :: left,right,total
  integer :: i,j,l
  do j=1,n
    do i=1,m
      total=0.0_dp
      do l=1,k
        if (transa=='n' .or. transa=='N') then
          left=a(i,l)
        else
          left=a(l,i)
        end if
        if (transb=='n' .or. transb=='N') then
          right=b(l,j)
        else
          right=b(j,l)
        end if
        total=total+left*right
      end do
      c(i,j)=alpha*total+beta*c(i,j)
    end do
  end do
end subroutine dgemm

module production_mrsf_extract
  implicit none
contains
"""
    return (
        preamble
        + "\n"
        + _extract_subroutine("mrsfcbc")
        + "\n"
        + _extract_subroutine("mrsfmntoia")
        + "\nend module production_mrsf_extract\n"
    )


@dataclass(frozen=True)
class ProductionResult:
    densities: dict[tuple[int, int], np.ndarray]
    sigma: dict[tuple[int, str], np.ndarray]


def _parse_fortran_output(output: str) -> ProductionResult:
    densities = {
        (multiplicity, column): np.zeros((7, NBF, NBF))
        for multiplicity in (1, 3)
        for column in range(len(PHYSICAL_LABELS[multiplicity]))
    }
    sigma = {
        (multiplicity, family): np.zeros(
            (NPACKED, len(PHYSICAL_LABELS[multiplicity]))
        )
        for multiplicity in (1, 3)
        for family in FAMILY_CASES
    }
    case_names = {value: name for name, value in FAMILY_CASES.items()}

    for line in output.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "D":
            multiplicity, column, channel, row, col = map(int, fields[1:6])
            densities[(multiplicity, column - 1)][
                channel - 1, row - 1, col - 1
            ] = float(fields[-1])
        elif fields[0] == "Y":
            multiplicity, case_id, row, column = map(int, fields[1:5])
            sigma[(multiplicity, case_names[case_id])][
                row - 1, column - 1
            ] = float(fields[-1])
    return ProductionResult(densities=densities, sigma=sigma)


@pytest.fixture(scope="module")
def production_result(tmp_path_factory: pytest.TempPathFactory) -> ProductionResult:
    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    if compiler is None:
        pytest.skip("GNU Fortran compiler is required for the production comparison")

    work = tmp_path_factory.mktemp("mrsf_spin_adapted_oracle")
    extracted = work / "production_mrsf_extract.F90"
    extracted.write_text(_production_extract_source())
    executable = work / "mrsf_spin_adapted_oracle"
    command = [
        compiler,
        "-std=f2018",
        "-O0",
        "-Wall",
        "-Wextra",
        "-J",
        str(work),
        str(extracted),
        str(CONVENTION_SOURCE),
        str(FORTRAN_DRIVER),
        "-o",
        str(executable),
    ]
    if platform.system() == "Darwin":
        command.extend(["-framework", "Accelerate"])
    else:
        command.append("-lblas")
    subprocess.run(command, cwd=work, check=True)
    output = subprocess.check_output([str(executable)], cwd=work, text=True)
    return _parse_fortran_output(output)


@pytest.mark.parametrize("multiplicity", [1, 3])
def test_explicit_seven_density_table_matches_production(
    production_result: ProductionResult,
    multiplicity: int,
) -> None:
    expected = _density_coefficients(multiplicity)
    actual = np.stack(
        [
            production_result.densities[(multiplicity, column)]
            for column in range(expected.shape[0])
        ]
    )
    np.testing.assert_allclose(actual, expected, atol=2.0e-14, rtol=0.0)


@pytest.mark.parametrize("multiplicity", [1, 3])
@pytest.mark.parametrize("family", list(FAMILY_CASES))
def test_full_signed_columns_match_independent_oracle(
    production_result: ProductionResult,
    multiplicity: int,
    family: str,
) -> None:
    _, expected = _oracle_sigma(multiplicity, family)
    actual = production_result.sigma[(multiplicity, family)]
    np.testing.assert_allclose(actual, expected, atol=5.0e-12, rtol=0.0)

    rows = _physical_rows(multiplicity)
    physical = actual[np.ix_(rows, range(actual.shape[1]))]
    np.testing.assert_allclose(
        physical, physical.T, atol=5.0e-12, rtol=0.0
    )


def test_oo_fold_and_triplet_forbidden_slots(
    production_result: ProductionResult,
) -> None:
    parent_resolved = np.array([[19.0, 17.0], [17.0, 19.0]])
    singlet = ISQRT2 * np.array([1.0, -1.0])
    triplet = ISQRT2 * np.array([1.0, 1.0])
    assert singlet @ parent_resolved @ singlet == pytest.approx(2.0)
    assert triplet @ parent_resolved @ triplet == pytest.approx(36.0)

    singlet_a0 = production_result.sigma[(1, "a0")]
    lr_column = PHYSICAL_LABELS[1].index("LR")
    assert abs(singlet_a0[PACKED_SLOTS["G"], lr_column]) > 1.0e-3
    assert abs(singlet_a0[PACKED_SLOTS["D"], lr_column]) > 1.0e-3

    for family in FAMILY_CASES:
        sigma = production_result.sigma[(3, family)]
        forbidden = [
            PACKED_SLOTS["G"],
            PACKED_SLOTS["D"],
            PACKED_SLOTS["R"],
        ]
        np.testing.assert_allclose(sigma[forbidden], 0.0, atol=1.0e-14)

    triplet_density = np.stack(
        [
            production_result.densities[(3, column)]
            for column in range(len(PHYSICAL_LABELS[3]))
        ]
    )
    assert np.max(np.abs(triplet_density[:, 6, O2, O1])) < 1.0e-14
    assert np.max(np.abs(triplet_density[:, 6, O1, O2])) < 1.0e-14


@pytest.mark.parametrize(
    ("family", "output_labels", "input_labels"),
    [
        ("coco", {"CO1", "CO2"}, {"CO1", "CO2"}),
        ("ovov", {"OV1", "OV2"}, {"OV1", "OV2"}),
        ("coov", {"CO1", "CO2", "OV1", "OV2"}, {"CO1", "CO2", "OV1", "OV2"}),
    ],
)
def test_spc_family_support_and_singlet_triplet_sign(
    production_result: ProductionResult,
    family: str,
    output_labels: set[str],
    input_labels: set[str],
) -> None:
    singlet = (
        production_result.sigma[(1, family)]
        - production_result.sigma[(1, "a0")]
    )
    triplet = (
        production_result.sigma[(3, family)]
        - production_result.sigma[(3, "a0")]
    )

    allowed_rows = {PACKED_SLOTS[label] for label in output_labels}
    allowed_s_columns = {
        index
        for index, label in enumerate(PHYSICAL_LABELS[1])
        if label in input_labels
    }
    for row in range(NPACKED):
        for column in range(singlet.shape[1]):
            if row not in allowed_rows or column not in allowed_s_columns:
                assert abs(singlet[row, column]) < 2.0e-14

    allowed_t_columns = {
        index
        for index, label in enumerate(PHYSICAL_LABELS[3])
        if label in input_labels
    }
    for row in range(NPACKED):
        for column in range(triplet.shape[1]):
            if row not in allowed_rows or column not in allowed_t_columns:
                assert abs(triplet[row, column]) < 2.0e-14

    common_columns = [
        PHYSICAL_LABELS[1].index(label)
        for label in PHYSICAL_LABELS[3][:5]
    ]
    np.testing.assert_allclose(
        triplet[:, :5], -singlet[:, common_columns], atol=5.0e-12, rtol=0.0
    )
    assert np.max(np.abs(singlet)) > 1.0e-3
    assert np.max(np.abs(triplet)) > 1.0e-3

    if family == "coov":
        co_rows = [PACKED_SLOTS["CO1"], PACKED_SLOTS["CO2"]]
        ov_rows = [PACKED_SLOTS["OV1"], PACKED_SLOTS["OV2"]]
        co_columns = [0, 1]
        ov_columns = [3, 4]
        assert np.max(np.abs(singlet[np.ix_(co_rows, ov_columns)])) > 1.0e-3
        assert np.max(np.abs(singlet[np.ix_(ov_rows, co_columns)])) > 1.0e-3
        np.testing.assert_allclose(
            singlet[np.ix_(co_rows, ov_columns)],
            singlet[np.ix_(ov_rows, co_columns)].T,
            atol=5.0e-12,
            rtol=0.0,
        )


@pytest.mark.parametrize(
    ("density_mutations", "output_mutations"),
    [
        (frozenset({"d5_plus"}), frozenset()),
        (frozenset({"d6_plus"}), frozenset()),
        (frozenset({"swap_channels_1_2"}), frozenset()),
        (frozenset({"swap_channels_3_4"}), frozenset()),
        (frozenset({"lr_input_wrong_sign"}), frozenset()),
        (frozenset({"lr_input_unnormalized"}), frozenset()),
        (frozenset(), frozenset({"d5_output_plus"})),
        (frozenset(), frozenset({"d6_output_plus"})),
        (frozenset(), frozenset({"lr_output_wrong_sign"})),
        (frozenset(), frozenset({"lr_output_unnormalized"})),
        (frozenset(), frozenset({"flip_channel_7"})),
        (frozenset(), frozenset({"retain_triplet_gd"})),
    ],
)
def test_negative_mutations_change_a_required_signed_column(
    density_mutations: frozenset[str],
    output_mutations: frozenset[str],
) -> None:
    largest_change = 0.0
    for multiplicity in (1, 3):
        for family in FAMILY_CASES:
            _, reference = _oracle_sigma(multiplicity, family)
            _, mutated = _oracle_sigma(
                multiplicity,
                family,
                density_mutations=density_mutations,
                output_mutations=output_mutations,
            )
            largest_change = max(
                largest_change, float(np.max(np.abs(reference - mutated)))
            )
    assert largest_change > 1.0e-3


def test_spc_channel_family_swap_is_detected() -> None:
    _, coco = _oracle_sigma(1, "coco")
    _, ovov = _oracle_sigma(1, "ovov")
    assert np.max(np.abs(coco - ovov)) > 1.0e-3
