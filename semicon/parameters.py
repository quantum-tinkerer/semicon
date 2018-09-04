import os
from types import SimpleNamespace

import yaml
import numpy as np
import pandas as pd
from scipy.constants import physical_constants as phys_const
from scipy.interpolate import interp1d


# General constants and globals
constants = {
    'm_0': phys_const['electron mass energy equivalent in MeV'][0] * 1e6,
    'phi_0': 2 * phys_const['mag. flux quantum'][0] * (1e9)**2,
    'mu_B': phys_const['Bohr magneton in eV/T'][0],
    'hbar': phys_const['Planck constant over 2 pi times c in MeV fm'][0]
}

taa = constants['hbar']**2 / 2 / constants['m_0']

# Load the parameters from the databank
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

_banks_names = ['winkler', 'lawaetz']

_parameters_names = [
    'E_0', 'Delta_0', 'P', 'm_c', 'gamma_1', 'gamma_2', 'gamma_3',
    'g_c',  'kappa', 'q'
]


def from_yaml(yml_str):
    pars = {k: v['parameters'] for k, v in yaml.load(yml_str).items()}
    df = pd.DataFrame(pars).T
    return df[_parameters_names]


def load_params(bankname):
    """Load material parameters from specified databank.

    output: pandas dataframe
    """
    if bankname not in _banks_names:
        msg = "Unkown bankname. Possible options are {}"
        raise ValueError(msg.format(_banks_names))
    fpath = os.path.join(BASE_DIR, 'databank', 'bank_' + bankname + '.yml')
    with open(fpath, 'r') as f:
        df = from_yaml(f.read())
    return df


# Renormalization of parameters
def renormalize_parameters(dict_pars, new_gamma_0=None,
                           bands=('gamma_6c', 'gamma_8v', 'gamma_7v')):
    """Renormalize parameters."""

    output = {}

    p = SimpleNamespace(**dict_pars)
    Ep = p.P**2 / taa
    output['P'] = p.P
    output['E_v'] = p.E_v
    output['E_0'] = p.E_0
    output['Delta_0'] = p.Delta_0

    # gamma_6c parameters
    if 'gamma_6c' in bands:
        # take special steps if the user want to renormalize gamma_0
        if new_gamma_0 is not None:
            if ('gamma_8v' not in bands) or ('gamma_7v' not in bands):
                raise ValueError('Cannot set different "gamma_0" '
                                 'without at least one hole band.')
            scale = 0
            if 'gamma_8v' in bands:
                scale += (2/3) / p.E_0

            if 'gamma_7v' in bands:
                scale += (1/3) / (p.E_0 + p.Delta_0)

            Ep = (p.gamma_0 - new_gamma_0) / scale
            output['P'] = np.sqrt(taa * Ep)
            output['gamma_0'] = new_gamma_0

        else:
            output['gamma_0'] = p.gamma_0
            if 'gamma_8v' in bands:
                output['gamma_0'] -= (2/3) * (Ep / p.E_0)

            if 'gamma_7v' in bands:
                output['gamma_0'] -= (1/3) * (Ep / (p.E_0 + p.Delta_0))

        # g-factor
        output['g_c'] = p.g_c
        if 'gamma_8v' in bands:
            output['g_c'] += (2/3) * (Ep / p.E_0)

        if 'gamma_7v' in bands:
            output['g_c'] -= (2/3) * (Ep / (p.E_0 + p.Delta_0))

    # now also for the gamma_7v and gamma_8v parameters
    if ('gamma_8v' in bands) or ('gamma_7v' in bands):
        output['gamma_1'] = p.gamma_1
        output['gamma_2'] = p.gamma_2
        output['gamma_3'] = p.gamma_3
        output['kappa'] = p.kappa
        output['q'] = p.q

        if 'gamma_6c' in bands:
            output['gamma_1'] -= (1/3) * (Ep / p.E_0)
            output['gamma_2'] -= (1/6) * (Ep / p.E_0)
            output['gamma_3'] -= (1/6) * (Ep / p.E_0)
            output['kappa'] -= (1/6) * (Ep / p.E_0)

    return output


# System specific parameter functions
def bulk(bank, material, new_gamma_0=None, new_P=None, valence_band_offset=0.0,
         bands=('gamma_6c', 'gamma_8v', 'gamma_7v'),
         extra_constants=None):
    """Get bulk parameters of a specified material."""

    if (new_gamma_0 is not None) and (new_P is not None):
        raise ValueError("'new_gamma_0' and 'new_P' cannot be used "
                         "simultaneously.")

    df_pars = load_params(bank)
    dict_pars = df_pars.loc[material].to_dict()
    dict_pars['gamma_0'] = 1 / dict_pars.pop('m_c')
    dict_pars['E_v'] = valence_band_offset

    if new_P is not None:
        dict_pars['P'] = new_P

    output = renormalize_parameters(dict_pars, new_gamma_0, bands)

    if extra_constants is not None:
        output.update(extra_constants)

    return output


def two_deg(parameters, widths, grid_spacing, extra_constants=None):
    """Get parameter functions for a specified 2D heterostructure.

    Parameters
    ----------
    parameters : sequence of dicts
        Material parameters for each material in the heterostructure.
        Only k.p parameters from each dictionary will be used.
    widths : sequence of numbers
        Width of each material in the heterostructure.
    grid_spacing : int, float
        Grid spacing that is used for discretization.
    extra_constants : dict
        Pass extra constants here.

    Returns
    -------
    parameters : dictionary of parameter functions
    walls : array of floats
    """

    def get_walls(a, Ws):
        walls = np.cumsum(Ws)[:-1] - 0.5 * a
        walls = np.insert(walls, 0, -a)
        walls = np.append(walls, sum(Ws))
        return walls

    def interp_sn_params(a, walls, values, parameter_name):
        xs = [x + d for x in walls[1:-1] for d in [-a/2, +a/2]]
        xs = [walls[0]] + xs + [walls[-1]]
        ys = [p[parameter_name] for p in values for i in range(2)]
        return interp1d(xs, ys, fill_value='extrapolate')

    # Varied parameters should probably be a union of available k·p parameters
    varied_parameters = ['E_0', 'E_v', 'Delta_0', 'P', 'kappa', 'g_c', 'q',
                         'gamma_0', 'gamma_1', 'gamma_2', 'gamma_3']

    walls = get_walls(grid_spacing, widths)

    output = {par: interp_sn_params(grid_spacing, walls, parameters, par)
              for par in varied_parameters}

    if extra_constants is not None:
        output.update(extra_constants)

    return output, walls


# Plotting helper function
def plot_2deg_bandedges(two_deg_params, xpos, walls=None, show_fig=False):
    """Plot band edges."""
    import matplotlib.pyplot as plt
    y1 = two_deg_params['E_v'](xpos)
    y2 = y1 + two_deg_params['E_0'](xpos)

    fig = plt.figure(figsize=(20, 5))
    plt.plot(xpos, y1, '-o')
    plt.plot(xpos, y2, '-o')

    if walls is not None:
        walls_y = [min([np.min(y1), np.min(y2)]),
                   max([np.max(y1), np.max(y2)])]
        for w in walls:
            plt.plot([w, w], walls_y, 'k--')

    if show_fig:
        plt.show()

    return fig
