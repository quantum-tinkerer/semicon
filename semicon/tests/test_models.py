import kwant
import sympy

from semicon.models.serialized import serialized_kp
from semicon.models.explicit_foreman import foreman as reference_foreman


def test_serialized_foreman():
    varied_parameters = ['E_0', 'E_v', 'Delta', 'P', 'kappa', 'g',
                     'gamma_0', 'gamma_1', 'gamma_2', 'gamma_3']

    substitutions = {v: v+'(x, y, z)' for v in varied_parameters}
    smp = kwant.continuum.sympify(serialized_kp, locals=substitutions)

    substitutions = {sympy.Symbol(k, commutative=False): kwant.continuum.sympify(v)
                 for k, v in substitutions.items()}
    foreman = reference_foreman.subs(substitutions)

    assert sympy.expand(smp - foreman) == sympy.zeros(8, 8)
