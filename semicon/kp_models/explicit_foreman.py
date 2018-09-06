"""
--------------------------------------------------------------------------
                Burt symmetrization of Kane Hamiltonian
--------------------------------------------------------------------------

This is explicit implementation of Hamiltonian in Sympy.
Semicon is using serialized version as a string, explicit implementation is
used in tests.


References:
[1] E. O. Kane, in Physics of III-V-Compounds, edited by R. K. Willardson and
    A. C. Beer, Vol. 1 of Semiconductors and Semimetals Academic Press,
    New York, London, 1666, p. 75.

[2] Bradley A. Foreman, Phys. Rev. B 56, R12748(R) - Published 15 November 1997
"""

import sympy
import sympy.physics
from sympy.physics.quantum import TensorProduct as kr

from .symbols import *

# Helper function definitions
def valence_term(i, j):
    kvec = [kx, ky, kz]
    """Return valence term in X, Y, Z basis."""
    if i == j:
        tmp = kvec[i] * L * kvec[i]
        for index in range(3):
            if index == i:
                continue
            else:
                tmp += kvec[index] * M * kvec[index]
        return tmp
    else:
        return kvec[i] * Np * kvec[j] + kvec[j] * Nm * kvec[i]


def spin_orbit():
    """Implementation of (4b) from PhysRevB 51, 16695."""
    data = sympy.zeros(8, 8)
    data[2, 1] = +sympy.I
    data[7, 1] = +1
    data[1, 2] = -sympy.I
    data[7, 2] = +sympy.I
    data[5, 3] = -1
    data[6, 3] = -sympy.I
    data[3, 5] = -1
    data[6, 5] = -sympy.I
    data[3, 6] = sympy.I
    data[5, 6] = sympy.I
    data[1, 7] = +1
    data[2, 7] = -sympy.I
    return (Delta / 3) * data


def get_foreman():
    # Basis transformation
    v_n = lambda i, N: sympy.Matrix([1 if i == n else 0 for n in range(N)])
    S, X, Y, Z = [v_n(i, 4) for i in range(4)]
    up, dw = [v_n(i, 2) for i in range(2)]

    X = sympy.I * X
    Y = sympy.I * Y
    Z = sympy.I * Z

    molenkamp_basis = [
        kr(up, S),
        kr(dw, S),
        +(1/sympy.sqrt(2)) * kr(up, X + sympy.I * Y),
        +(1/sympy.sqrt(6)) * (kr(dw, X + sympy.I * Y) - kr(up, 2*Z)),
        -(1/sympy.sqrt(6)) * (kr(up, X - sympy.I * Y) + kr(dw, 2*Z)),
        -(1/sympy.sqrt(2)) * kr(dw, X - sympy.I * Y),
        +(1/sympy.sqrt(3)) * (kr(dw, X + sympy.I * Y) + kr(up, Z)),
        +(1/sympy.sqrt(3)) * (kr(up, X - sympy.I * Y) - kr(dw, Z))
    ]

    subs_notation = {
        Ec: Ev + E0,
        L: -(g1 + 4 * g2) * (hbar**2/2/m0),
        M: -(g1 - 2 * g2) * (hbar**2/2/m0),
        Np: -(3 * g3 + (3 * kappa + 1)) * (hbar**2/2/m0),
        Nm: -(3 * g3 - (3 * kappa + 1)) * (hbar**2/2/m0),
        Ac: (g0 * hbar**2/2/m0)
    }

    Hs = spin_orbit()

    Hcc = sympy.Matrix([Ec + kx*Ac*kx + ky*Ac*ky + kz*Ac*kz])
    Hcv = +sympy.I * sympy.Matrix([[P*kx, P*ky, P*kz]])
    Hvc = -sympy.I * sympy.Matrix([kx*P, ky*P, kz*P])

    data = [[valence_term(i, j) for j in range(3)] for i in range(3)]
    Hvv = sympy.Matrix(data)

    H4 = sympy.BlockMatrix([[Hcc, Hcv], [Hvc, Hvv]])
    H8 = sympy.Matrix(sympy.BlockDiagMatrix(H4, H4)) + Hs

    dag = lambda x: x.conjugate().transpose()
    U_dag = sympy.Matrix(sympy.BlockMatrix([molenkamp_basis]))
    U = dag(U_dag)

    hamiltonian = (U * H8 * dag(U))
    hamiltonian = hamiltonian.subs(subs_notation)
    hamiltonian = (hamiltonian + (Ev - Delta / 3) *
                   sympy.diag(0, 0, 1, 1, 1, 1, 1, 1))
    return sympy.ImmutableMatrix(hamiltonian.expand())


foreman = get_foreman()
