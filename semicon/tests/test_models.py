import pytest

import numpy as np
import scipy.linalg as la
import kwant.continuum
import sympy


from semicon.misc import prettify, rotation_functionality_available
from semicon.models import Model
from semicon.kp_models import symbols


sigma_x = np.array(symbols.sigma_x.tolist(), dtype=complex)
sigma_y = np.array(symbols.sigma_y.tolist(), dtype=complex)
sigma_z = np.array(symbols.sigma_z.tolist(), dtype=complex)

Jx = np.array(symbols.Jx.tolist(), dtype=complex)
Jy = np.array(symbols.Jy.tolist(), dtype=complex)
Jz = np.array(symbols.Jz.tolist(), dtype=complex)


# Define helper functions
def isclose(a, b):
    if isinstance(a, sympy.matrices.MatrixBase):
        return prettify(a - b, zero_atol=1e-8) == sympy.zeros(*a.shape)
    else:
        return prettify(a - b, zero_atol=1e-8) == 0


# Baisc model tests

@pytest.mark.parametrize('ham_str', [
    "A_x * k_x**2 + A_y * k_y**2",
    "A_x * sigma_x * k_x**2 + A_y * sigma_y * k_y**2",
])
def test_model_creation(ham_str):
    smp = kwant.continuum.sympify(ham_str)
    assert smp == Model(ham_str).hamiltonian
    assert smp == Model(smp).hamiltonian

    locals = {'A_x': 1, 'A_y': 2}
    smp = kwant.continuum.sympify(ham_str, locals=locals)
    assert smp == Model(ham_str, locals=locals).hamiltonian


@pytest.mark.parametrize('str_input, str_output, decimals', [
    ("0.12345 * k_x**2", "0.123 * k_x**2", 3),
    ("0.12345 * k_x**2 * sigma_z", "0.123 * k_x**2 * sigma_z", 3),
    ("0.12345 * k_x**2", "0.1234 * k_x**2", 4),
    ("0.12345 * k_x**2 * sigma_z", "0.1234 * k_x**2 * sigma_z", 4),
])
def test_prettify_decimals(str_input, str_output, decimals):
    m = Model(str_input).prettify(decimals=decimals)
    assert isclose(m.hamiltonian, Model(str_output).hamiltonian)


@pytest.mark.parametrize('str_input, str_output, zero_atol', [
    ("1.123e-3 * k_x**2 + 5", "1.123e-3 * k_x**2 + 5", 1e-4),
    ("1.123e-3 * k_x**2 + 5", "5", 1e-2),
    ("(1.123e-3 * k_x**2 + 5) * sigma_z", "(1.123e-3 * k_x**2 + 5) * sigma_z", 1e-4),
    ("(1.123e-3 * k_x**2 + 5) * sigma_z", "5 * sigma_z", 1e-2),
])
def test_prettify_zero_atol(str_input, str_output, zero_atol):
    m = Model(str_input).prettify(zero_atol=zero_atol)
    assert isclose(m.hamiltonian, Model(str_output).hamiltonian)


@pytest.mark.parametrize('str_input, str_output, nsimplify', [
    ("1.7320508075688772 * k_x**2", "1.7320508075688772 * k_x**2", False),
    ("1.7320508075688772 * k_x**2", "sqrt(3) * k_x**2", True),
    ("1.7320508075688772 * k_x**2 + 5", "sqrt(3) * k_x**2 + 5", True),
    ("1.7320508075688772 * k_x**2 * sigma_z", "1.7320508075688772 * k_x**2 * sigma_z", False),
    ("1.7320508075688772 * k_x**2 * sigma_z", "sqrt(3) * k_x**2 * sigma_z", True),
])
def test_prettify_zero_nsimplify(str_input, str_output, nsimplify):
    m = Model(str_input).prettify(nsimplify=nsimplify)
    assert isclose(m.hamiltonian, Model(str_output).hamiltonian)


# Test rotation functionality

def test_spin_operators():
    model = Model("A_x * sigma_x * k_x**2 + A_y * sigma_y * k_y**2", spins=1/2)
    S = 0.5 * np.array([sigma_x, sigma_y, sigma_z])
    assert np.allclose(model.spin_operators, S)


R = np.array([[0, -1, 0],
              [+1, 0, 0],
              [0, 0, 1]])

def test_rotation():
    """Test rotation procedure.

    Here we test if general approach to implementing the Hamiltonian rotation
    is correct. It is we apply rotation explicitly to momenta and spins:
    1) assuming that k' = R.k and S' = R.S when Hamiltonian is defined using k and sigmas
    2) applying same to momenta but having spin include by working with matrices
    3) applying the builtin methods
    """
    if not rotation_functionality_available:
        return
    else:
        from scipy.spatial.transform import Rotation

    tmp_str = "alpha_{0} * k_{0} * sigma_{0}"
    ham_str = " + ".join(tmp_str.format(s) for s in ['x', 'y', 'z'])
    sx, sy, sz = sympy.symbols('sigma_x sigma_y sigma_z')


    ham1 = kwant.continuum.sympify(
        ham_str, locals={'sigma_x': sx, 'sigma_y': sy, 'sigma_z': sz}
    )

    ham2 = kwant.continuum.sympify(ham_str)

    S = 0.5 * np.array([sigma_x, sigma_y, sigma_z])

    get_subs = lambda R, v: {cprime: c for (cprime, c) in zip(v, R @ v)}
    subs = {
        **get_subs(R, sympy.Matrix(kwant.continuum.momentum_operators)),
        **get_subs(R, sympy.Matrix([sx, sy, sz]))
    }

    n = Rotation.from_dcm(R).as_rotvec()
    ns = np.sum([ni * si for (ni, si) in zip(n, S)], axis=0)
    U = la.expm(1j * ns)

    ham1_rotated = ham1.subs(subs, simultaneous=True).expand()

    ham2_substituted = ham2.subs(subs, simultaneous=True).expand()
    ham2_rotated = U @ ham2_substituted @ U.transpose().conjugate()

    a = kwant.continuum.sympify(str(ham1_rotated))
    b = kwant.continuum.sympify(str(ham2_rotated)).expand()

    assert isclose(a, b)
    assert isclose(a, Model(ham_str, spins=1/2).rotate(R).hamiltonian)
