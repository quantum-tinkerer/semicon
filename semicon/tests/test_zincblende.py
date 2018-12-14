import pytest


import kwant.continuum
import sympy

from semicon.models import ZincBlende
from semicon.misc import prettify

from semicon.kp_models.explicit_foreman import foreman as reference_foreman
from semicon.kp_models.explicit_zeeman import zeeman as reference_zeeman
from semicon.kp_models import symbols


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



# Sanity check of content: type, shape, included symbols...
# (if something is failing really badly, it should fail here)


@pytest.mark.parametrize('bands, shape', [
    ('gamma_6c', (2, 2)),
    ('gamma_8v', (4, 4)),
    ('gamma_7v', (2, 2)),
    (('gamma_6c', 'gamma_8v'), (6, 6)),
    (('gamma_6c', 'gamma_7v'), (4, 4)),
    (('gamma_8v', 'gamma_7v'), (6, 6)),
])
def test_hamiltonian_shape(bands, shape):
    model = ZincBlende(bands=bands)
    assert model.hamiltonian.shape == shape


@pytest.mark.parametrize('coords, components, must_have_symbols', [
    (None, 'foreman', kwant.continuum.momentum_operators),
    ('xyz', 'foreman', kwant.continuum.position_operators),
    (None, 'zeeman', symbols.magnetic_symbols),
])
def test_hamiltonian_symbols(coords, components, must_have_symbols):
    model = ZincBlende(components=components, parameter_coords=coords)
    atoms = model.hamiltonian.atoms(sympy.Symbol)
    assert set(must_have_symbols).issubset(atoms)
