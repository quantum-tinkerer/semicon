import pandas as pd
import os

from scipy.constants import physical_constants
from scipy.interpolate import interp1d
from types import SimpleNamespace
import numpy as np

import matplotlib.pyplot as plt


######  general constants and globals
c = physical_constants['speed of light in vacuum'][0]
val_hbar = physical_constants['Planck constant over 2 pi in eV s'][0]
val_m0 = physical_constants['electron mass energy equivalent in MeV'][0]
val_m0 = val_m0 / (c*10**9)**2 * 10**6
val_mu_B = physical_constants['Bohr magneton in eV/T'][0]
val_phi_0 = physical_constants['mag. flux quantum'][0] * (10**9)**2
taa = val_hbar**2 / 2.0 / val_m0

constants = {
    'm_0': val_m0,
    'phi_0': val_phi_0,
    'mu_B': val_mu_B,
    'hbar': val_hbar
}

###### loading parameters from databank
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

_banks_names = ['winkler', 'lawaetz']
def load_params(bankname):
    """Load material parameters from specified databank.

    output: pandas dataframe
    """
    if bankname not in _banks_names:
        msg = "Unkown bankname. Possible options are {}"
        raise ValueError(msg.format(_banks_names))
    fname = 'bank_' + bankname + '.csv'
    fpath = os.path.join(BASE_DIR, 'databank', fname)
    return pd.read_csv(fpath, index_col=0)


###### renormalization of parameters
def renormalize_parameters(dict_pars, new_gamma_0=None):
    """Renormalize parameters"""

    param = SimpleNamespace(**dict_pars)

    e0_plus_d0 = param.E_0 + param.Delta_0
    e0_plus_23d0 = param.E_0 + (2 / 3) * param.Delta_0

    scale1 = param.E_0 * e0_plus_d0 / e0_plus_23d0

    if new_gamma_0 is not None:
        scale1 = taa * (param.gamma_0 - new_gamma_0) * scale1
        param.P = np.sqrt(scale1)
        param.gamma_0 = new_gamma_0
    else:
        param.gamma_0 -= (param.P**2 / taa) / scale1

    scale2 = (param.P**2 / taa) / param.E_0

    param.gamma_1 -= (1 / 3) * scale2
    param.gamma_2 -= (1 / 6) * scale2
    param.gamma_3 -= (1 / 6) * scale2
    param.kappa -= (1 / 6) * scale2

    param.g_c += (2 / 3) * scale2 * (param.Delta_0 / e0_plus_d0)

    return param.__dict__


###### system specific parameter functions
def bulk(bank, material, new_gamma_0=None):
    """Get bulk parameters of a specified material."""
    df_pars = load_params(bank)
    dict_pars = df_pars.loc[material].to_dict()
    dict_pars['gamma_0'] = 1 / dict_pars.pop('m_c')
    return renormalize_parameters(dict_pars, new_gamma_0)


def two_deg(bank, materials, widths, valence_band_offsets, grid_spacing,
            new_gamma_0=None, extra_constants=None):
    """Get parameter functions for a specified 2D system.


    To do:
    - way provide parameters, not only data "bank" name,
    - specify which parameters should be varied (default all defined parameters)
    """

    def get_walls(a, Ws):
        walls = np.array([sum(Ws[:(i+1)]) for i in range(len(Ws)-1)]) - 0.5 * a
        walls = np.insert(walls, 0, -a)
        walls = np.append(walls, sum(Ws))
        return walls

    def interp_sn_params(a, walls, values, parameter_name):
        xs = [x+d for x in walls[1:-1] for d in [-a/2, +a/2]]
        xs = [walls[0]] + xs + [walls[-1]]
        ys = [p[parameter_name] for p in values for i in range(2)]
        return interp1d(xs, ys, fill_value='extrapolate')

    # varied parameters should probably be a union of available k.p parameters,
    varied_parameters = ['E_0', 'E_v', 'Delta_0', 'P', 'kappa', 'g_c',
                         'gamma_0', 'gamma_1', 'gamma_2', 'gamma_3']

    df_pars = load_params(bank)

    parameters = []
    for mat_name, offset in zip(materials, valence_band_offsets):
        dict_pars = df_pars.loc[mat_name].to_dict()
        dict_pars['gamma_0'] = 1 / dict_pars.pop('m_c')
        dict_pars['E_v'] = offset
        dict_pars = renormalize_parameters(dict_pars, new_gamma_0)
        parameters.append(dict_pars)

    walls = get_walls(grid_spacing, widths)

    output = {}
    for par_name in varied_parameters:
        output[par_name] = interp_sn_params(grid_spacing, walls, parameters, par_name)

    if extra_constants is not None:
        output.update(extra_constants)

    return output, walls


###### helper plotting function
def plot_2deg_bandedges(two_deg_params, xpos, walls=None, show_fig=False):
    """Plot band edges """
    y1 = two_deg_params['E_v'](xpos)
    y2 = y1 + two_deg_params['E_0'](xpos)

    fig = plt.figure(figsize=(20, 5))
    plt.plot(xpos, y1, '-o')
    plt.plot(xpos, y2, '-o')

    if walls is not None:
        walls_y = [min([np.min(y1), np.min(y2)]), max([np.max(y1), np.max(y2)])]
        for w in walls:
            plt.plot([w, w], walls_y, 'k--')

    if show_fig:
        plt.show()

    return fig
