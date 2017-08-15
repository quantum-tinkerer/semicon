import os
import json
import kwant

import sympy


# parameters varied in k.p Hamiltonian
varied_parameters = ['E_0', 'E_v', 'Delta_0', 'P', 'kappa', 'g',
                     'gamma_0', 'gamma_1', 'gamma_2', 'gamma_3',
                     'g', 'mu_B', 'kappa', 'q']


##### read cache
def _load_cache():
    """Load cached models.

    File semicon/kp_models/cache.json should be created on package build.
    """
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    fname = os.path.join(BASE_DIR, 'kp_models', 'cache.json')
    with open(fname) as f:
        models_cache = json.load(f)
    return models_cache

_models_cache = _load_cache()



##### module functions
def validate_coords(coords):
    """Validate coords in the same way it happens in kwant.continuum."""
    coords = list(coords)
    if coords != sorted(coords):
        raise ValueError("The argument 'coords' must be sorted.")
    if any(c not in 'xyz' for c in coords):
        raise ValueError("The argument 'coords' may only contain "
                         "'x', 'y', or 'z'.")
    return coords


def foreman(coords=None, components=('foreman',)):
    """Return 8x8 k.p Hamiltonian following Burt-Foreman symmetrization.

    Parameters
    ----------
    coords : sequence of strings
        Spatial dependents of parameters, e.g. ``coords='xyz'``
    components : sequence of strings
        k.p components, .e.g. ``components=['foreman', 'zeeman']``

    Returns
    -------
    kp_hamiltonian : sympy object
    """
    if coords is not None:
        coords = validate_coords(coords)
        str_coords ='({})'.format(", ".join(coords))
        subs = {v: v + str_coords for v in varied_parameters}
    else:
        subs = {}

    hamiltonian_components = [
        kwant.continuum.sympify(_models_cache[c], locals=subs)
        for c in components
    ]

    return sympy.ImmutableMatrix(sympy.MatAdd(*hamiltonian_components))
