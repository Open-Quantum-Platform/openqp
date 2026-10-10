"""Matrix manipulation utilities"""
from math import sqrt, floor
import numpy as np


def unpack_symmetric(matrix_packed, matrix_dim=None):
    """Convert symmetric packed matrix to full square format"""
    # Guess dimension:
    if matrix_dim:
        dim = matrix_dim
    else:
        # Compute 'triangular root' of matrix size
        tr_dim = len(matrix_packed)
        dim = floor((sqrt(8 * tr_dim + 1) - 1) / 2)
    # Repack to numpy 2D array
    matrix_square = np.zeros([dim, dim])
    k = 0
    for i in range(dim):
        matrix_square[i, 0:i + 1] = matrix_packed[k:k + i + 1]
        k += i + 1
    matrix_square += np.tril(matrix_square, -1).T
    return matrix_square


def pack_matrix(matrix):
    """
       Saving the lower triangular part of the matrix in a packed format,
       which corresponds to the upper triangular part in Fortran notation.
    """
    n, m = matrix.shape
    return matrix[np.tril_indices(n)]


def orb_to_dens(v, x):
    """
      Compute density matrix from a set of orbitals and respective occupation numbers.
      Compute the transformation: D = V * diag(X) * V^T
      Args:
       V  - matrix of orbitals
       X  - vector of occupation numbers
      Returns:
       D  - density matrix
    """
    occup = len(x)
    # Compute D = V * diag(X) * V^T
    d = np.einsum('ki,k,kj -> ij', v[:occup, :], x, v[:occup, :])

    return pack_matrix(d)

class DispersionModel:
    # empty class to skip dftd4 calculation

    def __init__(self, atoms, coordinates):
        self.natom = len(atoms)

    def get_dispersion(self, model, grad):
        corr = {
            'energy': 0,
            'gradient': np.zeros((self.natom, 3))
        }

        return corr


def DampingParam(method=None):
    # empty function to skip dftd4 calculation
    return None
