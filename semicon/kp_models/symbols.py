import sympy
import sympy.physics
from sympy import Matrix, sqrt

# Momenta and positions
momentum_symbols = kx, ky, kz = sympy.symbols("k_x k_y k_z", commutative=False)
position_symbols = sympy.symbols("x y z", commutative=False)


# Symbols
Ec, Ac, P, M, L = sympy.symbols("E_c A_c P M L", commutative=False)
N, Np, Nm = sympy.symbols("N N_+ N_-", commutative=False)

# Gamma parameters
E0, Ev, g0 = sympy.symbols("E_0 E_v gamma_0", commutative=False)
g1, g2, g3 = sympy.symbols("gamma_1 gamma_2 gamma_3", commutative=False)
Delta, kappa = sympy.symbols("Delta_0 kappa", commutative=False)
P, mu, gbar = sympy.symbols("P mu gammabar", commutative=False)
hbar, m0 = sympy.symbols("hbar, m_0")


# ************** symbols **************
mu_b, g, kappa, q = sympy.symbols("mu_B, g_c, kappa, q", commutative=False)
magnetic_symbols = Bx, By, Bz = sympy.symbols("B_x, B_y, B_z")


# ************** Magnetic field Hamiltonian (Winkler form) **************

# Turn off black formatting to visually compare the matrices to a text book source
# fmt: off
Tx = 1 / (3 * sqrt(2)) * Matrix([
    [-sqrt(3), 0, 1, 0],
    [0, -1, 0, sqrt(3)]
])

Ty = -sympy.I / (3 * sqrt(2)) * Matrix([
    [sqrt(3), 0, 1, 0],
    [0, 1, 0, sqrt(3)]
])

Tz = sqrt(2) / 3 * Matrix([
    [0, 1, 0, 0],
    [0, 0, 1, 0]
])

Txx = 1 / (3 * sqrt(2)) * Matrix([
    [0, -1, 0, sqrt(3)],
    [-sqrt(3), 0, 1, 0]
])

Tyy = 1 / (3 * sqrt(2)) * Matrix([
    [0, -1, 0, -sqrt(3)],
    [sqrt(3), 0, 1, 0]
])

Tzz = sqrt(2) / 3 * Matrix([
    [0, 1, 0, 0],
    [0, 0, -1, 0]
])

Tyz = sympy.I / (2 * sqrt(6)) * Matrix([
    [-1, 0, -sqrt(3), 0],
    [0, sqrt(3), 0, 1]
])

Tzx = 1 / (2 * sqrt(6)) * Matrix([
    [-1, 0, sqrt(3), 0],
    [0, sqrt(3), 0, -1]
])

Txy = sympy.I / sqrt(6) * Matrix([
    [0, 0, 0, -1],
    [-1, 0, 0, 0]
])

Jx = sympy.Rational(1, 2) * Matrix([
    [0, sqrt(3), 0, 0],
    [sqrt(3), 0, 2, 0],
    [0, 2, 0, sqrt(3)],
    [0, 0, sqrt(3), 0]
])

Jy = (
    sympy.I
    * sympy.Rational(1, 2)
    * Matrix([
        [0, -sqrt(3), 0, 0],
        [sqrt(3), 0, -2, 0],
        [0, 2, 0, -sqrt(3)],
        [0, 0, sqrt(3), 0],
    ])
)
# fmt: on

Jz = sympy.Rational(1, 2) * sympy.diag(3, 1, -1, -3)

sigma_0 = sympy.eye(2)
sigma_x = sympy.physics.matrices.msigma(1)
sigma_y = sympy.physics.matrices.msigma(2)
sigma_z = sympy.physics.matrices.msigma(3)
