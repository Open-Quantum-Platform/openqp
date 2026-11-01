"""Tag record element counts must be computed in 64-bit integers.

tagarray v1.0.0 multiplied a record's extents in a 32-bit ``int``.  A record
with more than 2^31 - 1 elements, such as the dense nbf^4 AO ERI tensor of
``int2e`` at nbf >= 216, got a negative or wrapped count, so tagarray
allocated no memory (or too little) but still reported success; ``int2e``
then crashed on its first write.  ``external/fix_tagarray_count64.py`` carries
the upstream fix until a tagarray release contains it.

The probe shape (2^32 + 3, 2^61 + 1) has more elements than int64_t can
count, so a correct library refuses it before allocating anything
(DATA_INCORRECT_DIMENSIONS) and ``oqp_alloc`` returns -3 (create refused).
A real allocation failure (MEMORY_ALLOCATION_ERROR) maps to -2 instead; see
``test_native_create_status_maps_to_alloc_code`` in
``test_oqpdata_setitem_errors.py``.  A 32-bit
accumulator wraps it to 3 elements, so a library with the defect allocates 24
bytes and reports success.  Neither case allocates real memory, so the test is
safe on every CI runner.
"""
import numpy as np
import pytest

TA_DIMENSIONS_LENGTH = 12


@pytest.fixture
def data():
    from oqp.molecule.oqpdata import OQPData
    return OQPData()


def _alloc(data, key, shape):
    from oqp import ffi, lib
    type_id = data.typemap.index(np.float64)
    req = ffi.new("char []", key.encode("ascii"))
    tid = ffi.new("int32_t *", type_id)
    dims = ffi.new(f"int64_t[{TA_DIMENSIONS_LENGTH}]", list(shape))
    ndims = ffi.new("int32_t *", len(shape))
    ptr = ffi.new("void **")
    return lib.oqp_alloc(data._data, req, tid, ndims, dims, ptr), ptr[0]


def test_count_beyond_32_bits_is_not_wrapped(data):
    from oqp import ffi
    shape = ((1 << 32) + 3, (1 << 61) + 1)
    n, ptr = _alloc(data, "OQP::PROBE_COUNT64", shape)
    # 3 is what a 32-bit accumulator produces for this shape.
    assert n != 3, "tagarray wrapped the element count to 32 bits"
    # oqp_alloc returns -3 when tagarray refuses the create.
    assert n == -3 and ptr == ffi.NULL, (n, ptr)


def test_small_asymmetric_record_round_trips(data):
    value = np.arange(105, dtype=np.float64).reshape(3, 5, 7) * 0.5 - 7.0
    data["OQP::PROBE_COUNT64_SMALL"] = value
    back = np.asarray(data["OQP::PROBE_COUNT64_SMALL"])
    assert back.shape == (3, 5, 7)
    assert np.array_equal(back, value)
