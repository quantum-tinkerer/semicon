import sympy
import sympy.physics.matrices

from .symbols import *


def kill_complex(smp):
    subs = {s.conjugate(): s for s in smp.atoms(sympy.Symbol)}
    return smp.subs(subs)


B_6c6c = sympy.Rational(1, 2) * g * mu_b * (sigma_x * Bx +
                                           sigma_y * By + sigma_z * Bz)
B_8v8v = -2 * mu_b * (kappa * (Jx * Bx + Jy * By + Jz * Bz) +
                     q * (Jx**3 * Bx + Jy**3 * By + Jz**3 * Bz))
B_7v8v = -3 * mu_b * kappa * (Tx * Bx + Ty * By + Tz * Bz)
B_8v7v = -3 * mu_b * kappa * \
    (Tx * Bx + Ty * By + Tz * Bz).conjugate().transpose()
B_7v7v = -2 * kappa * mu_b * (sigma_x * Bx + sigma_y * By + sigma_z * Bz)


zeeman = sympy.BlockMatrix([[B_6c6c, sympy.zeros(2, 4), sympy.zeros(2, 2)],
                            [sympy.zeros(4, 2), B_8v8v, B_8v7v],
                            [sympy.zeros(2, 2), B_7v8v, B_7v7v]]).as_explicit()


# adjust symbols that got conjugated
zeeman = kill_complex(zeeman)
