import kwant
import sympy

from semicon.models import foreman
from semicon.kp_models.explicit_foreman import foreman as reference_foreman


def test_serialized_foreman():
    smp = foreman('xyz')

    varied_parameters = ['E_0', 'E_v', 'Delta_0', 'P', 'kappa', 'g',
                     'gamma_0', 'gamma_1', 'gamma_2', 'gamma_3']

    substitutions = {v: v+'(x, y, z)' for v in varied_parameters}
    substitutions = {sympy.Symbol(k, commutative=False): kwant.continuum.sympify(v)
                     for k, v in substitutions.items()}
    reference = reference_foreman.subs(substitutions)

    assert sympy.expand(smp - reference) == sympy.zeros(8, 8)
