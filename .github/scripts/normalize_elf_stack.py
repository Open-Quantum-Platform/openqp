#!/usr/bin/env python3
"""Remove executable-stack requests from packaged OpenQP ELF libraries."""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import io
import struct
import sysconfig
from pathlib import Path


PT_GNU_STACK = 0x6474E551
PF_X = 0x1


def _elf_layout(data: bytes | bytearray) -> tuple[str, int, int, int, int]:
    if len(data) < 16 or data[:4] != b"\x7fELF":
        raise ValueError("not an ELF file")

    elf_class = data[4]
    byte_order = data[5]
    if byte_order == 1:
        endian = "<"
    elif byte_order == 2:
        endian = ">"
    else:
        raise ValueError(f"unsupported ELF byte order: {byte_order}")

    if elf_class == 1:
        if len(data) < 52:
            raise ValueError("truncated ELF32 header")
        phoff = struct.unpack_from(endian + "I", data, 28)[0]
        phentsize = struct.unpack_from(endian + "H", data, 42)[0]
        phnum = struct.unpack_from(endian + "H", data, 44)[0]
        flags_offset = 24
        minimum_entry_size = 28
    elif elf_class == 2:
        if len(data) < 64:
            raise ValueError("truncated ELF64 header")
        phoff = struct.unpack_from(endian + "Q", data, 32)[0]
        phentsize = struct.unpack_from(endian + "H", data, 54)[0]
        phnum = struct.unpack_from(endian + "H", data, 56)[0]
        flags_offset = 4
        minimum_entry_size = 8
    else:
        raise ValueError(f"unsupported ELF class: {elf_class}")

    if phnum == 0xFFFF:
        raise ValueError("extended ELF program-header counts are unsupported")
    if phentsize < minimum_entry_size:
        raise ValueError(f"invalid ELF program-header size: {phentsize}")
    if phoff + phentsize * phnum > len(data):
        raise ValueError("ELF program-header table extends past end of file")
    return endian, phoff, phentsize, phnum, flags_offset


def gnu_stack_flags(data: bytes | bytearray) -> int:
    endian, phoff, phentsize, phnum, flags_offset = _elf_layout(data)
    for index in range(phnum):
        entry = phoff + index * phentsize
        segment_type = struct.unpack_from(endian + "I", data, entry)[0]
        if segment_type == PT_GNU_STACK:
            return struct.unpack_from(endian + "I", data, entry + flags_offset)[0]
    raise ValueError("ELF file has no PT_GNU_STACK program header")


def clear_executable_stack(path: Path) -> bool:
    data = bytearray(path.read_bytes())
    endian, phoff, phentsize, phnum, flags_offset = _elf_layout(data)
    for index in range(phnum):
        entry = phoff + index * phentsize
        segment_type = struct.unpack_from(endian + "I", data, entry)[0]
        if segment_type != PT_GNU_STACK:
            continue
        offset = entry + flags_offset
        flags = struct.unpack_from(endian + "I", data, offset)[0]
        if flags & PF_X:
            struct.pack_into(endian + "I", data, offset, flags & ~PF_X)
            path.write_bytes(data)
            if gnu_stack_flags(path.read_bytes()) & PF_X:
                raise RuntimeError(f"failed to clear executable stack: {path}")
            return True
        return False
    raise ValueError(f"ELF file has no PT_GNU_STACK program header: {path}")


def openqp_package_libraries() -> list[Path]:
    package_lib = Path(sysconfig.get_path("purelib")) / "oqp" / "lib"
    if not package_lib.is_dir():
        raise FileNotFoundError(f"installed OpenQP library directory missing: {package_lib}")
    libraries = sorted(
        path for path in package_lib.glob("*.so*")
        if path.is_file() and path.read_bytes()[:4] == b"\x7fELF"
    )
    if not libraries:
        raise FileNotFoundError(f"no installed OpenQP ELF libraries found: {package_lib}")
    return libraries


def _record_hash(path: Path) -> tuple[str, str]:
    payload = path.read_bytes()
    encoded = base64.urlsafe_b64encode(hashlib.sha256(payload).digest())
    return f"sha256={encoded.rstrip(b'=').decode('ascii')}", str(len(payload))


def update_openqp_record(paths: list[Path], purelib: Path | None = None) -> None:
    """Update and verify installed wheel RECORD entries for normalized files."""
    package_root = (purelib or Path(sysconfig.get_path("purelib"))).resolve()
    records = sorted(
        child / "RECORD"
        for child in package_root.glob("*.dist-info")
        if child.name.lower().startswith("openqp-") and (child / "RECORD").is_file()
    )
    if len(records) != 1:
        raise RuntimeError(f"expected one installed OpenQP RECORD, found {records}")
    record = records[0]
    rows = list(csv.reader(record.read_text(encoding="utf-8").splitlines()))
    row_map = {row[0]: row for row in rows if len(row) == 3}

    expected: dict[str, tuple[str, str]] = {}
    for path in paths:
        resolved = path.resolve()
        try:
            relative = resolved.relative_to(package_root).as_posix()
        except ValueError as exc:
            raise RuntimeError(f"normalized file is outside site-packages: {path}") from exc
        if relative not in row_map:
            raise RuntimeError(f"OpenQP RECORD has no entry for {relative}")
        expected[relative] = _record_hash(resolved)
        row_map[relative][1], row_map[relative][2] = expected[relative]

    output = io.StringIO(newline="")
    csv.writer(output, lineterminator="\n").writerows(rows)
    record.write_text(output.getvalue(), encoding="utf-8")

    verified_rows = {
        row[0]: (row[1], row[2])
        for row in csv.reader(record.read_text(encoding="utf-8").splitlines())
        if len(row) == 3
    }
    for relative, expected_record in expected.items():
        if verified_rows.get(relative) != expected_record:
            raise RuntimeError(f"failed to update OpenQP RECORD entry: {relative}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    parser.add_argument(
        "--openqp-package",
        action="store_true",
        help="normalize every ELF shared library in the installed oqp/lib directory",
    )
    parser.add_argument(
        "--update-record",
        action="store_true",
        help="update and verify installed OpenQP dist-info/RECORD entries",
    )
    args = parser.parse_args()
    paths = list(args.paths)
    if args.openqp_package:
        paths.extend(openqp_package_libraries())
    if not paths:
        parser.error("provide a library path or --openqp-package")

    seen = set()
    normalized_paths = []
    for path in paths:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        normalized_paths.append(resolved)
        changed = clear_executable_stack(resolved)
        state = "cleared" if changed else "already non-executable"
        print(f"PT_GNU_STACK {state}: {resolved}")
    if args.update_record:
        update_openqp_record(normalized_paths)
        print("Updated and verified OpenQP dist-info/RECORD")


if __name__ == "__main__":
    main()
