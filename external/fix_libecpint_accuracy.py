#!/usr/bin/env python3
"""Make the libecpint 1.0.7 integrals accurate enough for analytic gradients.

Two defects each leave the ECP integrals with relative errors of 1e-6 to 1e-2.
Because the integral and its nuclear derivative are evaluated from different
angular-momentum shells, the errors do not cancel: the analytic ECP gradient
of HOBr/LANL2DZ disagrees with finite differences by 1e-5 Ha/bohr.

1. The code generator writes the unrolled type-2 angular coefficients with the
   default stream precision, so the compiled Q*.cpp literals keep only six
   significant digits (16*pi^2 = 157.9136704... becomes 157.914).
2. The one-point Gauss-Chebyshev rule stops at the first level whose estimate
   passes its convergence test.  Two coarse levels whose abscissae all miss a
   narrow integrand pass by coincidence, and the type-1 radial integrals keep
   errors of order 1e-2.  Every caller tabulates the integrand on all
   abscissae beforehand, so summing the full grid costs only the additions;
   the test still sets the convergence flag the callers act on.
"""
from pathlib import Path
import re
import sys

# Upstream lines carry trailing blanks, so match each line up to optional
# trailing whitespace instead of byte-exactly.
GEN_INC = (re.compile(r"#include <algorithm>[ \t]*\n#include <set>"),
           "#include <algorithm>\n#include <iomanip>\n#include <limits>\n#include <set>")
GEN_PREC = (re.compile(r"(\tstd::ofstream outfile\(ofname\);[ \t]*\n)"),
            r"\1\t// OpenQP: write the angular coefficients at full double precision.\n"
            r"\toutfile << std::setprecision(std::numeric_limits<double>::max_digits10);\n")
GEN_DONE = "std::setprecision(std::numeric_limits<double>::max_digits10)"

QUAD = (re.compile(
    r"while \(n < maxN && !converged\) \{(?P<body>[ \t]*\n"
    r"[ \t]*// Compute T_\{2n\+1\}[ \t]*\n"
    r"[ \t]*T2n1 = Tn \+ sumTerms\(f, params, n, start, end, p, 2\);[ \t]*\n"
    r"[ \t]*\n"
    r"[ \t]*// Check convergence[ \t]*\n"
    r"[ \t]*dT = T2n1 - 2\.0\*Tn;[ \t]*\n"
    r"[ \t]*n = 2\*n \+ 1;[ \t]*\n)"
    r"(?P<ind>[ \t]*)if \(dT\*dT <= fabs\(T2n1 - Tn12\)\*tolerance\) \{[ \t]*\n"
    r"[ \t]*converged = true;[ \t]*\n"
    r"[ \t]*\} else \{[ \t]*\n"
    r"[ \t]*Tn12 = 4\.0 \* Tn;[ \t]*\n"
    r"[ \t]*Tn = T2n1;[ \t]*\n"
    r"[ \t]*p /= 2;[ \t]*\n"
    r"[ \t]*\}[ \t]*\n"),
    "while (n < maxN) { // OpenQP: always sum the full (tabulated) grid"
    "\\g<body>"
    "\\g<ind>if (dT*dT <= fabs(T2n1 - Tn12)*tolerance) converged = true;\n"
    "\\g<ind>Tn12 = 4.0 * Tn;\n"
    "\\g<ind>Tn = T2n1;\n"
    "\\g<ind>p /= 2;\n")
QUAD_DONE = "while (n < maxN) { // OpenQP: always sum the full (tabulated) grid"


def apply(path, edits, done, label):
    text = path.read_text()
    if text.count(done) == 1:
        print(f"[OpenQP] libecpint {label}: already applied")
        return
    for pattern, repl in edits:
        text, n = pattern.subn(repl, text)
        if n != 1:
            raise RuntimeError(f"Unrecognized libecpint source for {label}: {path}")
    path.write_text(text)
    print(f"[OpenQP] libecpint {label}: applied")


def patch(root):
    root = Path(root)
    apply(root / "src" / "generate.cpp", [GEN_INC, GEN_PREC], GEN_DONE,
          "full-precision generated angular coefficients")
    apply(root / "src" / "lib" / "gaussquad.cpp", [QUAD], QUAD_DONE,
          "full-grid one-point Gauss-Chebyshev quadrature")


if __name__ == "__main__":
    patch(sys.argv[1])
