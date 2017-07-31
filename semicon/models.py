from .kp_models.serialized import serialized_kp
import kwant


# parameters varied in k.p Hamiltonian
varied_parameters = ['E_0', 'E_v', 'Delta_0', 'P', 'kappa', 'g',
                     'gamma_0', 'gamma_1', 'gamma_2', 'gamma_3']


def validate_coords(coords):
    """Validate coords in the same way it happens in kwant.continuum."""
    coords = list(coords)
    if coords != sorted(coords):
        raise ValueError("The argument 'coords' must be sorted.")
    if any(c not in 'xyz' for c in coords):
        raise ValueError("The argument 'coords' may only contain "
                         "'x', 'y', or 'z'.")
    return coords


def kane(coords):
    coords = validate_coords(coords)
    str_coords ='({})'.format(", ".join(coords))

    subs = {v: v + str_coords for v in varied_parameters}
    return kwant.continuum.sympify(serialized_kp, locals=subs)
