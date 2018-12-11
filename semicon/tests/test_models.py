import pytest

import kwant
import sympy

from semicon.misc import prettify
from semicon.models import Model, ZincBlende
from kp_models.explicit_foreman import foreman as reference_foreman
from kp_models.explicit_zeeman import zeeman as reference_zeeman



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
