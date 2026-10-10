#!/usr/bin/env python3
"""Embed an RTC2 tables .bin into a C++ TU so the library is self-contained.
Usage: embed_tables.py data/routec_tables.bin src/routec_tables_embedded.cpp"""
import sys
src, dst = sys.argv[1], sys.argv[2]
data = open(src, "rb").read()
assert data[:4] == b"RTC2", data[:4]
n = len(data)
with open(dst, "w") as f:
    f.write(f"// AUTO-GENERATED from {src} (RTC2) -- DO NOT EDIT.\n")
    f.write(f"// Regenerate: python3 tools/embed_tables.py {src} {dst}\n")
    f.write("// Embedding the rotation tables makes libopenqp_gpu.so self-contained: no\n")
    f.write("// external routec_tables.bin needed at runtime (OQP_ROUTEC_TABLES overrides).\n")
    f.write("#include <cstddef>\n")
    f.write(f"extern const unsigned char routec_tables_embedded[{n}];\n")
    f.write("extern const size_t routec_tables_embedded_len;\n")
    f.write(f"const size_t routec_tables_embedded_len = {n};\n")
    f.write(f"const unsigned char routec_tables_embedded[{n}] = {{\n")
    for i in range(0, n, 20):
        f.write("".join(f"{b:d}," for b in data[i:i+20]) + "\n")
    f.write("};\n")
print(f"embedded {n} bytes -> {dst}")
