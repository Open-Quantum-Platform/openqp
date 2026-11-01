#!/usr/bin/env python3
"""Match the v1.0.0 record deleter to its 64-byte aligned allocation."""
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


if __name__ == "__main__":
    patch_alignment(sys.argv[1])
