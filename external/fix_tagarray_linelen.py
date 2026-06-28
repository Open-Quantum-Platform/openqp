#!/usr/bin/env python3
# OpenQP build patch: tagarray v1.0.0 ships container.f90 with a 144-char
# `use, intrinsic :: iso_c_binding, only: ...` line. tagarray's own CMake
# compiles the library with -std=f2018, whose 132-char free-form limit
# truncates that line -> c_null_ptr/c_associated drop out -> cascade of
# "has no IMPLICIT type" / "not a member" errors under gfortran >= 12.
# This wraps the single offending line with continuations (<=132 chars,
# fully f2018-conformant). Idempotent: a no-op once wrapped or if upstream
# fixes it. Mirror this fix upstream in Open-Quantum-Platform/tagarray.
import sys, pathlib

src = pathlib.Path(sys.argv[1]) / "source" / "API" / "Fortran" / "container.f90"
if not src.is_file():
    print(f"[OpenQP] tagarray patch: {src} not found, skipping")
    sys.exit(0)

text = src.read_text()
long_line = ("  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, "
             "c_float, c_double, c_double_complex, c_ptr, c_f_pointer, "
             "c_null_ptr, c_associated")
wrapped = ("  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_float, &\n"
           "                                         c_double, c_double_complex, c_ptr, &\n"
           "                                         c_f_pointer, c_null_ptr, c_associated")

if long_line in text:
    src.write_text(text.replace(long_line, wrapped))
    print("[OpenQP] tagarray patch: wrapped container.f90 iso_c_binding line")
else:
    print("[OpenQP] tagarray patch: already conformant, no change")
