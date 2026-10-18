"""Shared container plumbing for the characteristic-net stages (MATLAB value semantics)."""

from __future__ import annotations

import numpy as np

from ._matlab import GrowArray, GrowColumns, GrowVector


class Net:
    """Mutable working copy of (X, Y, U, V, LENG_INDI, patch indices).

    Each stage function copies its inputs into a ``Net`` (MATLAB passes arrays
    by value), works with 1-based indices exactly like the MATLAB code, and
    returns NumPy copies.
    """

    def __init__(self, X, Y, U, V, LENG_INDI=None, i_tri=None, j_tri=None, i_quad=None, j_quad=None):
        self.X = GrowArray(X)
        self.Y = GrowArray(Y)
        self.U = GrowArray(U)
        self.V = GrowArray(V)
        self.LENG_INDI = GrowVector(LENG_INDI)
        self.i_tri = GrowColumns(3, i_tri)
        self.j_tri = GrowColumns(3, j_tri)
        self.i_quad = GrowColumns(4, i_quad)
        self.j_quad = GrowColumns(4, j_quad)

    def point(self, J, I):
        return self.X[J, I], self.Y[J, I], self.U[J, I], self.V[J, I]

    def mirror(self, J, I):
        """Point (J, I) reflected about the axis (y -> -y, v -> -v)."""
        return self.X[J, I], -self.Y[J, I], self.U[J, I], -self.V[J, I]

    def set_point(self, J, I, x, y, u, v):
        self.X[J, I] = x
        self.Y[J, I] = y
        self.U[J, I] = u
        self.V[J, I] = v

    def outputs(self, indI, ind_tri, ind_quad):
        return (
            indI,
            self.X.to_numpy(),
            self.Y.to_numpy(),
            self.U.to_numpy(),
            self.V.to_numpy(),
            self.LENG_INDI.to_numpy(),
            ind_tri,
            ind_quad,
            self.i_tri.to_numpy(),
            self.j_tri.to_numpy(),
            self.i_quad.to_numpy(),
            self.j_quad.to_numpy(),
        )


def as_int(value) -> int:
    return int(np.real(value))
