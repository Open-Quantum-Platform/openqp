"""openqp-gpu: Python interface to the GPU density-fitting SCF library.

    from openqp_gpu import solve_rhf, build_df_tensor
"""
from .driver import solve_rhf, solve_uhf
from .build_df import build_df_tensor

__all__ = ["solve_rhf", "solve_uhf", "build_df_tensor"]
__version__ = "0.1.0"
