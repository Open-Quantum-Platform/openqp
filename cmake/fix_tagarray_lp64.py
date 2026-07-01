#!/usr/bin/env python3
"""Add LP64-friendly reserve_data overloads to tagarray v0.0.6."""

from pathlib import Path
import sys


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: fix_tagarray_lp64.py <tagarray-source-dir>", file=sys.stderr)
        return 2

    container = Path(sys.argv[1]) / "source" / "API" / "Fortran" / "container.f90"
    text = container.read_text()
    if "reserve_data_default" in text:
        return 0

    text = text.replace(
        "    procedure, public  :: reserve_data\n",
        "    generic, public :: reserve_data => reserve_data_i64, reserve_data_default\n"
        "    procedure, private :: reserve_data_i64\n"
        "    procedure, private :: reserve_data_default\n",
    )
    text = text.replace(
        "  subroutine reserve_data(this, tag, datatype, array_size, array_shape, options, comment)\n",
        "  subroutine reserve_data_i64(this, tag, datatype, array_size, array_shape, options, comment)\n",
    )
    text = text.replace(
        "  end subroutine reserve_data\n",
        """  end subroutine reserve_data_i64
  subroutine reserve_data_default(this, tag, datatype, array_size, array_shape, options, comment)
    class(container_t), intent(inout) :: this
    character(kind=TA_CHAR, len=*),           intent(in) :: tag
    integer(c_int32_t),                             intent(in) :: datatype
    integer,                                        intent(in) :: array_size
    integer,                              optional, intent(in) :: array_shape(:)
    integer,                              optional, intent(in) :: options(TA_OPTIONS_LENGTH)
    character(kind=TA_CHAR, len=*), optional, intent(in) :: comment
    !
    integer(c_int64_t), allocatable :: array_shape_64(:)
    integer(c_int64_t) :: options_64(TA_OPTIONS_LENGTH)

    if (present(options)) then
      options_64 = int(options, c_int64_t)
    endif

    if (present(array_shape)) then
      array_shape_64 = int(array_shape, c_int64_t)
      if (present(options)) then
        call reserve_data_i64(this, tag, datatype, int(array_size, c_int64_t), array_shape_64, options_64, comment)
      else
        call reserve_data_i64(this, tag, datatype, int(array_size, c_int64_t), array_shape_64, comment=comment)
      endif
    else
      if (present(options)) then
        call reserve_data_i64(this, tag, datatype, int(array_size, c_int64_t), options=options_64, comment=comment)
      else
        call reserve_data_i64(this, tag, datatype, int(array_size, c_int64_t), comment=comment)
      endif
    endif
  end subroutine reserve_data_default
""",
    )
    container.write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
