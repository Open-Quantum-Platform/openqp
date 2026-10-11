#!/usr/bin/env python3
"""Build-consumed TagArray 1.0 allocation and ownership corrections."""
from pathlib import Path
import sys


def patch_alignment(root):
    header = Path(root) / "include" / "tagarray" / "hidden" / "Record.hpp"
    text = header.read_text()
    old = "operator delete[](ptr, std::align_val_t(32));"
    new = "operator delete[](ptr, std::align_val_t(64));"
    if text.count(old) == 1 and text.count(new) == 1:
        header.write_text(text.replace(old, new))
    elif old in text or text.count(new) != 2:
        raise RuntimeError("Unrecognized TagArray record allocation/deallocation")
    print("[OpenQP] TagArray record allocation and deletion both use alignment 64")


def patch_creation(root):
    # Container owns only fully allocated records. A null payload is not success.
    header = Path(root) / "include" / "tagarray" / "Container.hpp"
    text = header.read_text()
    old = """    this->records_.insert(std::make_pair(
        tag, new hidden::Record(type_id, dimensions, data, description)));
    return defines::OK;"""
    new = """    try {
      auto record = std::make_shared<hidden::Record>(type_id, dimensions, data,
                                                      description);
      if (record->byte_count() > 0 && record->raw_data<uint8_t *>() == nullptr)
        return defines::MEMORY_ALLOCATION_ERROR;
      this->records_.insert(std::make_pair(tag, std::move(record)));
    } catch (const std::bad_alloc &) {
      return defines::MEMORY_ALLOCATION_ERROR;
    }
    return defines::OK;"""
    if text.count(old) == 1 and new not in text:
        text = text.replace(old, new)
        text = text.replace("#include <memory>", "#include <memory>\n#include <new>")
        header.write_text(text)
    elif old in text or text.count(new) != 1 or "#include <new>" not in text:
        raise RuntimeError("Unrecognized TagArray record creation")
    print("[OpenQP] TagArray rejects failed payload allocation before insertion")


if __name__ == "__main__":
    patch_alignment(sys.argv[1])
    patch_creation(sys.argv[1])
