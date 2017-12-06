import sympy
import sympy.physics.matrices


# ************** symbols **************
g, kappa, q = sympy.symbols("g_c, kappa, q", commutative=False)
u_b, Bx, By, Bz = sympy.symbols('mu_B, B_x, B_y, B_z', real=True)


# ************** Magnetic field Hamiltonian (Winkler form) **************
Tx = 1 / (3 * sympy.sqrt(2)) * sympy.Matrix([[-sympy.sqrt(3), 0, 1, 0],
                                             [0, -1, 0, sympy.sqrt(3)]])

Ty = -sympy.I / (3 * sympy.sqrt(2)) * sympy.Matrix([[sympy.sqrt(3), 0, 1, 0],
                                                    [0, 1, 0, sympy.sqrt(3)]])

Tz = sympy.sqrt(2) / 3 * sympy.Matrix([[0, 1, 0, 0], [0, 0, 1, 0]])


Txx = 1 / (3 * sympy.sqrt(2)) * sympy.Matrix([[0, -1, 0, sympy.sqrt(3)],
                                              [-sympy.sqrt(3), 0, 1, 0]])

Tyy = 1 / (3 * sympy.sqrt(2)) * sympy.Matrix([[0, -1, 0, -sympy.sqrt(3)],
                                              [sympy.sqrt(3), 0, 1, 0]])

Tzz = sympy.sqrt(2) / 3 * sympy.Matrix([[0, 1, 0, 0],
                                        [0, 0, -1, 0]])

Tyz = sympy.I / (2 * sympy.sqrt(6)) * sympy.Matrix([[-1, 0, -sympy.sqrt(3), 0],
                                                    [0, sympy.sqrt(3), 0, 1]])

Tzx = 1 / (2 * sympy.sqrt(6)) * sympy.Matrix([[-1, 0, sympy.sqrt(3), 0],
                                              [0, sympy.sqrt(3), 0, -1]])

Txy = sympy.I / sympy.sqrt(6) * sympy.Matrix([[0, 0, 0, -1],
                                                [-1, 0, 0, 0]])

sigma_0 = sympy.eye(2)
sigma_x = sympy.physics.matrices.msigma(1)
sigma_y = sympy.physics.matrices.msigma(2)
sigma_z = sympy.physics.matrices.msigma(3)

Jx = sympy.Rational(1, 2) * sympy.Matrix([[0, sympy.sqrt(3), 0, 0],
                                          [sympy.sqrt(3), 0, 2, 0],
                                          [0, 2, 0, sympy.sqrt(3)],
                                          [0, 0, sympy.sqrt(3), 0]])

Jy = sympy.I * sympy.Rational(1, 2) * sympy.Matrix([[0, -sympy.sqrt(3), 0, 0],
                                                    [sympy.sqrt(3), 0, -2, 0],
                                                    [0, 2, 0, -sympy.sqrt(3)],
                                                    [0, 0, sympy.sqrt(3), 0]])

Jz = sympy.Rational(1, 2) * sympy.diag(3, 1, -1, -3)


B_6c6c = sympy.Rational(1, 2) * g * u_b * (sigma_x * Bx +
                                           sigma_y * By + sigma_z * Bz)
B_8v8v = -2 * u_b * (kappa * (Jx * Bx + Jy * By + Jz * Bz) +
                     q * (Jx**3 * Bx + Jy**3 * By + Jz**3 * Bz))
B_7v8v = -3 * u_b * kappa * (Tx * Bx + Ty * By + Tz * Bz)
B_8v7v = -3 * u_b * kappa * \
    (Tx * Bx + Ty * By + Tz * Bz).conjugate().transpose()
B_7v7v = -2 * kappa * u_b * (sigma_x * Bx + sigma_y * By + sigma_z * Bz)


zeeman = sympy.BlockMatrix([[B_6c6c, sympy.zeros(2, 4), sympy.zeros(2, 2)],
                            [sympy.zeros(4, 2), B_8v8v, B_8v7v],
                            [sympy.zeros(2, 2), B_7v8v, B_7v7v]]).as_explicit()


u_b, Bx, By, Bz = sympy.symbols('mu_B B_x B_y B_z')
zeeman_subs = {sympy.Symbol(s.name, real=True): s for s in [u_b, Bx, By, Bz]}
zeeman = zeeman.subs(zeeman_subs)
