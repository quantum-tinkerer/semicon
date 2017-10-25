import kwant
import sympy

from semicon.models import foreman
from semicon.kp_models.explicit_foreman import foreman as reference_foreman
from semicon.kp_models.explicit_zeeman import zeeman as reference_zeeman


def test_serialized_foreman():
    smp = foreman('xyz')

    varied_parameters = ['E_0', 'E_v', 'Delta_0', 'P', 'kappa', 'g_c', 'q',
                         'gamma_0', 'gamma_1', 'gamma_2', 'gamma_3']

    substitutions = {v: v+'(x, y, z)' for v in varied_parameters}
    substitutions = {sympy.Symbol(k, commutative=False): kwant.continuum.sympify(v)
                     for k, v in substitutions.items()}
    reference = reference_foreman.subs(substitutions)

    assert sympy.expand(smp - reference) == sympy.zeros(8, 8)


def test_serialized_zeeman():
    smp = foreman('xyz', components=['zeeman'])

    varied_parameters = ['E_0', 'E_v', 'Delta_0', 'P', 'kappa', 'g_c', 'q',
                         'gamma_0', 'gamma_1', 'gamma_2', 'gamma_3']

    substitutions = {v: v+'(x, y, z)' for v in varied_parameters}
    substitutions = {sympy.Symbol(k, commutative=False): kwant.continuum.sympify(v)
                     for k, v in substitutions.items()}
    reference = reference_zeeman.subs(substitutions)

    assert sympy.expand(smp - reference) == sympy.zeros(8, 8)


def test_serialized_foremanzeeman():
    smp = foreman('xyz', components=['zeeman', 'foreman'])

    varied_parameters = ['E_0', 'E_v', 'Delta_0', 'P', 'kappa', 'g_c', 'q',
                         'gamma_0', 'gamma_1', 'gamma_2', 'gamma_3']

    substitutions = {v: v+'(x, y, z)' for v in varied_parameters}
    substitutions = {sympy.Symbol(k, commutative=False): kwant.continuum.sympify(v)
                     for k, v in substitutions.items()}
    reference = (
        reference_foreman.subs(substitutions)
        + reference_zeeman.subs(substitutions)
    )

    assert sympy.expand(smp - reference) == sympy.zeros(8, 8)
