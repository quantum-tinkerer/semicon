import kwant
import sympy

from semicon.misc import prettify
from semicon.models import ZincBlende
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


# Test functions
def test_serialized_foreman():
    smp = ZincBlende(parameter_coords='xyz')
    assert isclose(smp.hamiltonian, reference_foreman)


def test_serialized_zeeman():
    smp = ZincBlende(components=['zeeman'], parameter_coords='xyz')
    assert isclose(smp.hamiltonian, reference_zeeman)


def test_serialized_foremanzeeman():
    smp = ZincBlende(components=['zeeman', 'foreman'], parameter_coords='xyz')
    assert isclose(smp.hamiltonian, reference_foreman + reference_zeeman)
