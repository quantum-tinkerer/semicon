# -*- coding: utf-8 -*-
import os
import json
import kwant

import sympy

from .kp_models.symbols import sigma_x, sigma_y, sigma_z, Jx, Jy, Jz


# Parameters varied in the k·p Hamiltonian
varied_parameters = ['E_0', 'E_v', 'Delta_0', 'P', 'kappa', 'g_c', 'q',
                     'gamma_0', 'gamma_1', 'gamma_2', 'gamma_3']


# Read the cache
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


# Module functions

def validate_coords(coords):
    """Validate coords in the same way it happens in kwant.continuum."""
    coords = list(coords)
    if coords != sorted(coords):
        raise ValueError("The argument 'coords' must be sorted.")
    if any(c not in 'xyz' for c in coords):
        raise ValueError("The argument 'coords' may only contain "
                         "'x', 'y', or 'z'.")
    return coords


def _band_indices(bands):
    band_indices = {
        'gamma_6c': [0, 1],
        'gamma_8v': [2, 3, 4, 5],
        'gamma_7v': [6, 7]
    }

    for b in bands:
        if b not in band_indices:
            raise ValueError("{} is not a proper band".format(b))

    indices = []
    for band in bands:
        indices += band_indices[band]
    return indices


def foreman(coords=None, components=('foreman',),
            bands=('gamma_6c', 'gamma_8v', 'gamma_7v')):
    """Return 8x8 k·p Hamiltonian following Burt-Foreman symmetrization.

    Parameters
    ----------
    coords : sequence of strings
        Spatial dependents of parameters, e.g. ``coords='xyz'``
    components : sequence of strings
        k·p components, e.g. ``components=['foreman', 'zeeman']``
    bands : sequence of strings
        k·p bands, e.g. ``bands=['gamma_6c']

    Returns
    -------
    kp_hamiltonian : sympy object
    """
    if coords is not None:
        coords = validate_coords(coords)
        str_coords = '({})'.format(", ".join(coords))
        subs = {v: v + str_coords for v in varied_parameters}
    else:
        subs = {}

    hamiltonian_components = [
        kwant.continuum.sympify(_models_cache[c], locals=subs)
        for c in components
    ]

    hamiltonian = sympy.ImmutableMatrix(sympy.MatAdd(*hamiltonian_components))

    # if the default bands are selected, we just return the "hamiltonian"
    if tuple(bands) == ('gamma_6c', 'gamma_8v', 'gamma_7v'):
        return hamiltonian

    indices = _band_indices(bands)
    return hamiltonian[:, indices][indices, :]


def spin_operators(bands=('gamma_6c', 'gamma_8v', 'gamma_7v')):
    Sx = sympy.BlockDiagMatrix(sigma_x / 2, Jx, sigma_x / 2).as_explicit()
    Sy = sympy.BlockDiagMatrix(sigma_y / 2, Jy, sigma_y / 2).as_explicit()
    Sz = sympy.BlockDiagMatrix(sigma_z / 2, Jz, sigma_z / 2).as_explicit()

    indices = _band_indices(bands)
    output = [S[:, indices][indices, :] for S in [Sx, Sy, Sz]]
    return output
