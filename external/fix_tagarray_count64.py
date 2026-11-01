#!/usr/bin/env python3
"""Compute tagarray v1.0.0 record element counts in int64_t.

tagarray v1.0.0 evaluates a record's element count as
``std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<int64_t>())``.
The initial value ``1`` is an ``int``, so the running product is stored back in
a 32-bit ``int`` after every multiplication.  A record with more than 2^31 - 1
elements (the dense nbf^4 AO ERI tensor at nbf >= 216) then gets a negative or
wrapped count: no memory, or too little memory, is allocated, yet
``Container::create`` returns OK and the caller writes through a null or
undersized data pointer.

This applies the upstream fix (Open-Quantum-Platform/tagarray, branch
fix/int64-record-count): the count is computed in int64_t with an overflow
check, and ``Container::create`` reports DATA_INCORRECT_DIMENSIONS or
MEMORY_ALLOCATION_ERROR instead of storing an unallocated record.  Remove this
script when _OQP_TAGARRAY_VERSION moves to a release that contains that fix.
"""
from pathlib import Path
import sys

RECORD_EDITS = [
    ("#include <cstdint>\n", "#include <cstdint>\n#include <limits>\n"),
    ("""public:
  Record(const int32_t type_id,""", """public:
  // Element count of a record with extents `dims`, evaluated in int64_t.
  // Returns -1 for a negative extent, or when the element count or the byte
  // count (count * itemsize) does not fit in int64_t.
  inline static int64_t checked_count(const std::vector<int64_t> &dims,
                                      const int32_t itemsize) noexcept {
    int64_t count = 1;
    for (const auto extent : dims) {
      if (extent < 0)
        return -1;
      if (extent != 0 && count > std::numeric_limits<int64_t>::max() / extent)
        return -1;
      count *= extent;
    }
    if (itemsize > 0 && count > std::numeric_limits<int64_t>::max() / itemsize)
      return -1;
    return count;
  }

  Record(const int32_t type_id,"""),
    ("""    this->count_ =
        std::accumulate(std::begin(this->dims_), std::end(this->dims_), 1,
                        std::multiplies<int64_t>());
""", """    this->count_ = checked_count(this->dims_, this->itemsize_);
"""),
    ("""    assert(this->count_ == std::accumulate(std::begin(shape), std::end(shape),
                                           1, std::multiplies<int64_t>()));
""", """    assert(this->count_ == checked_count(shape, this->itemsize_));
"""),
    ("""    int64_t shift = std::accumulate(std::begin(shape), std::end(shape), 1,
                                    std::multiplies<int64_t>());
    int64_t new_count = shift + this->count_;
""", """    int64_t shift = checked_count(shape, this->itemsize_);
    if (shift < 0 ||
        shift > std::numeric_limits<int64_t>::max() / this->itemsize_ -
                    this->count_)
      return defines::DATA_INCORRECT_DIMENSIONS;
    int64_t new_count = shift + this->count_;
"""),
]

CONTAINER_EDITS = [
    ("""    this->records_.insert(std::make_pair(
        tag, new hidden::Record(type_id, dimensions, data, description)));
    return defines::OK;""", """    auto record = std::make_shared<hidden::Record>(type_id, dimensions, data,
                                                   description);
    if (record->count() < 0)
      return defines::DATA_INCORRECT_DIMENSIONS;
    if (record->byte_count() > 0 && record->raw_data<uint8_t *>() == nullptr)
      return defines::MEMORY_ALLOCATION_ERROR;
    this->records_.insert(std::make_pair(tag, record));
    return defines::OK;"""),
]


def apply(path, edits, marker):
    text = path.read_text()
    if marker in text:
        return False
    for old, new in edits:
        if text.count(old) != 1:
            raise RuntimeError(f"Unrecognized TagArray source in {path}")
        text = text.replace(old, new)
    path.write_text(text)
    return True


def patch_count64(root):
    include = Path(root) / "include" / "tagarray"
    changed = apply(include / "hidden" / "Record.hpp", RECORD_EDITS,
                    "checked_count(")
    changed |= apply(include / "Container.hpp", CONTAINER_EDITS,
                     "return defines::MEMORY_ALLOCATION_ERROR;")
    if "std::accumulate" in (include / "hidden" / "Record.hpp").read_text():
        raise RuntimeError("TagArray record count still uses std::accumulate")
    state = "patched" if changed else "already patched"
    print(f"[OpenQP] TagArray record counts use int64_t ({state})")


if __name__ == "__main__":
    patch_count64(sys.argv[1])
