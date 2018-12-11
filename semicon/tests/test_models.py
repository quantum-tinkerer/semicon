import pytest

import numpy as np
import scipy.linalg as la
from scipy.spatial.transform import Rotation

import kwant
import sympy


from semicon.misc import prettify
from semicon.models import Model, ZincBlende

from kp_models.explicit_foreman import foreman as reference_foreman
from kp_models.explicit_zeeman import zeeman as reference_zeeman
from kp_models.symbols import sigma_x, sigma_y, sigma_z, Jx, Jy, Jz


sigma_x = np.array(sigma_x.tolist(), dtype=complex)
sigma_y = np.array(sigma_y.tolist(), dtype=complex)
sigma_z = np.array(sigma_z.tolist(), dtype=complex)

Jx = np.array(Jx.tolist(), dtype=complex)
Jy = np.array(Jy.tolist(), dtype=complex)
Jz = np.array(Jz.tolist(), dtype=complex)


# Prepare reference Hamiltonian with proper commutivities
varied_parameters = ['E_0', 'E_v', 'Delta_0', 'P', 'kappa', 'g_c', 'q',
                     'gamma_0', 'gamma_1', 'gamma_2', 'gamma_3']

substitutions = {v: v+'(x, y, z)' for v in varied_parameters}

reference_foreman = kwant.continuum.sympify(
    str(reference_foreman),
    locals=substitutions
)

reference_zeeman = kwant.continuum.sympify(
    str(reference_zeeman),
    locals=substitutions
)


# Define helper functions
def isclose(a, b):
    return prettify(a - b, zero_atol=1e-8) == sympy.zeros(*a.shape)


# Test correctness of serialization process
# main point here is to check if order of operators is preserved as it is
# important in context of proper symmetrization.
def test_serialized_foreman():
    smp = ZincBlende(parameter_coords='xyz')
    assert isclose(smp.hamiltonian, reference_foreman)


def test_serialized_zeeman():
    smp = ZincBlende(components=['zeeman'], parameter_coords='xyz')
    assert isclose(smp.hamiltonian, reference_zeeman)


def test_serialized_foremanzeeman():
    smp = ZincBlende(components=['zeeman', 'foreman'], parameter_coords='xyz')
    assert isclose(smp.hamiltonian, reference_foreman + reference_zeeman)



# Sanity check of returned types (if something breaks just a little, it will
# break here for sure...)

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
