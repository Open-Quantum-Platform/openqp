"""A tag write that cannot be stored must raise, never print and return.

``OQPData.__setitem__`` used to wrap the dtype lookup and the native
``oqp_alloc`` call in a bare ``except:`` that printed "The type of your data
is not allowed" and returned.  Every failed tag write was then silent, so a
caller such as the SCF stability safeguard's state restore continued with the
old contents of the tag.  These tests require that an unsupported dtype raises
TypeError naming the tag and dtype, and that a failing ``oqp_alloc``
propagates as an exception, while supported writes still round-trip.
"""
import numpy as np
import pytest


@pytest.fixture
def data():
    from oqp.molecule.oqpdata import OQPData
    return OQPData()


@pytest.mark.parametrize("value", [
    np.array([True, False]),
    np.array(["a", "bc"]),
    np.array([1.0, 2.0], dtype=np.float16),
    np.array([None, 1], dtype=object),
])
def test_unsupported_dtype_raises_type_error(data, value):
    key = "OQP::PROBE_UNSUPPORTED_DTYPE"
    with pytest.raises(TypeError) as err:
        data[key] = value
    message = str(err.value)
    assert key in message
    assert value.dtype.str in message


def test_supported_dtypes_round_trip(data):
    # Non-symmetric values so a transposed or reordered copy is visible.
    f64 = np.arange(6, dtype=np.float64).reshape(2, 3) * 1.5 - 2.0
    i64 = np.array([3, -1, 7], dtype=np.int64)
    i32 = np.array([[1, 2], [3, 4]], dtype=np.int32)
    data["OQP::PROBE_F64"] = f64
    data["OQP::PROBE_I64"] = i64
    data["OQP::PROBE_I32"] = i32
    data["OQP::PROBE_STR"] = "probe"
    np.testing.assert_array_equal(np.asarray(data["OQP::PROBE_F64"]), f64)
    np.testing.assert_array_equal(np.asarray(data["OQP::PROBE_I64"]), i64)
    np.testing.assert_array_equal(np.asarray(data["OQP::PROBE_I32"]), i32)
    assert data["OQP::PROBE_STR"] == "probe"


class _FailingAlloc:
    """Delegate to the real cffi lib but make ``oqp_alloc`` report failure."""

    def __init__(self, real, status):
        self._real = real
        self._status = status

    def oqp_alloc(self, *args):
        return self._status

    def __getattr__(self, name):
        return getattr(self._real, name)


@pytest.mark.parametrize("status, exc", [
    (-2, MemoryError), (-1, RuntimeError), (-3, ValueError)])
def test_failed_native_alloc_propagates(data, monkeypatch, status, exc):
    import oqp.molecule.oqpdata as oqpdata_mod

    key = "OQP::PROBE_ALLOC_FAILURE"
    data[key] = np.array([1.0, 2.0])
    monkeypatch.setattr(oqpdata_mod, "lib", _FailingAlloc(oqpdata_mod.lib, status))
    with pytest.raises(exc) as err:
        data[key] = np.array([5.0, 6.0, 7.0])
    assert key in str(err.value)
    monkeypatch.undo()
    # The failed write left the previous contents in place; the caller saw
    # an exception instead of silently continuing with them.
    np.testing.assert_array_equal(np.asarray(data[key]), [1.0, 2.0])


def _tag_status(data, key):
    """Raw ``oqp_get`` return for ``key``: >= 0 when a tag record exists."""
    from oqp import ffi, lib
    from oqp.molecule.oqpdata import TA_DIMENSIONS_LENGTH

    return lib.oqp_get(data._data, ffi.new('char []', bytes(key, 'ascii')),
                       ffi.new('int32_t *'), ffi.new('int32_t *'),
                       ffi.new(f'int64_t[{TA_DIMENSIONS_LENGTH}]'),
                       ffi.new('void **'))


@pytest.mark.parametrize("key, value, structs", [
    ("tda", True, ("tddft",)),
    ("debug_mode", True, ("tddft", "mpiinfo")),
    ("usempi", False, ("mpiinfo",)),
    ("charge", -1, ("mol_prop",)),
])
def test_struct_member_write_is_not_stored_as_tag(data, key, value, structs):
    # A struct member (a bool here, which tagarray cannot store) is written to
    # every native struct that has it and the call returns; it must neither
    # raise after the member changed nor leave a stray tag of the same name.
    data[key] = value
    for name in structs:
        assert getattr(getattr(data._data, name), key) == value
    assert data[key] == value
    assert _tag_status(data, key) < 0


def test_refused_native_create_raises_before_copy(data):
    # tagarray rejects a tag containing a control character.  oqp_alloc must
    # report that refusal itself (-3) rather than look the tag up afterwards.
    key = "OQP::PROBE\x01REFUSED"
    with pytest.raises(ValueError) as err:
        data[key] = np.array([1.0, 2.0, 3.0])
    assert "refused" in str(err.value)
    assert _tag_status(data, key) < 0


# tagarray v1.0.0 status values (defines.csv of the pinned release).  The
# harness takes the statuses from the tagarray_defines module itself; these
# values only make a renumbering in a tagarray upgrade visible here.
_TA_STATUS_CODES = [
    ("TA_OK", 0, 0),
    ("TA_MEMORY_ALLOCATION_ERROR", -1, -2),
    ("TA_DATA_INCORRECT_DIMENSIONS", -10, -3),
    ("TA_INCORRECT_TAG", -14, -3),
    ("TA_CONTAINER_RECORD_EXISTS", -4, -3),
    ("TA_DATA_TOO_MANY_DIMENSIONS", -8, -3),
    ("TA_UNKNOWN_ERROR", -1000, -3),
    ("TA_NOT_IMPLEMENTED", 1000, -3),
]


def test_native_create_status_maps_to_alloc_code():
    # oqp_alloc must report a tagarray allocation failure as -2 (MemoryError
    # in OQPData) and keep -3 (ValueError) for an invalid tag, invalid
    # dimensions or a record conflict.  A real allocation failure cannot be
    # induced portably (GCC 15 libstdc++ on macOS traps in a failing nothrow
    # aligned new[]), so the Fortran harness drives the mapping that
    # oqp_alloc uses with the named tagarray status constants.
    from oqp import ffi, lib

    nmax = 16
    statuses = ffi.new(f"int32_t[{nmax}]")
    codes = ffi.new(f"int64_t[{nmax}]")
    n = ffi.new("int64_t *")
    lib.oqp_alloc_status_selftest(nmax, statuses, codes, n)
    assert n[0] == len(_TA_STATUS_CODES)
    got = [(statuses[i], codes[i]) for i in range(n[0])]
    want = [(status, code) for _, status, code in _TA_STATUS_CODES]
    assert got == want, list(zip((name for name, _, _ in _TA_STATUS_CODES), got))

    # Too small an output buffer is refused, not overrun.
    lib.oqp_alloc_status_selftest(2, statuses, codes, n)
    assert n[0] == -1


@pytest.mark.parametrize("_name, status, code", _TA_STATUS_CODES[1:])
def test_alloc_code_reaches_python_as_documented_exception(
        data, monkeypatch, _name, status, code):
    # Second half of the route: the code oqp_alloc returns for each failing
    # tagarray status raises the matching OQPData exception.
    import oqp.molecule.oqpdata as oqpdata_mod

    exc = MemoryError if code == -2 else ValueError
    monkeypatch.setattr(oqpdata_mod, "lib", _FailingAlloc(oqpdata_mod.lib, code))
    with pytest.raises(exc):
        data["OQP::PROBE_STATUS_ROUTE"] = np.array([1.0, 2.0])
