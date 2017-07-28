import pandas as pd
import os

from scipy.constants import physical_constants
from types import SimpleNamespace
import numpy as np



######  general constants and globals
c = physical_constants['speed of light in vacuum'][0]
val_hbar = physical_constants['Planck constant over 2 pi in eV s'][0]
val_m0 = physical_constants['electron mass energy equivalent in MeV'][0]
val_m0 = val_m0 / (c*10**9)**2 * 10**6
val_mu_B = physical_constants['Bohr magneton in eV/T'][0]
val_phi_0 = physical_constants['mag. flux quantum'][0] * (10**9)**2
taa = val_hbar**2 / 2.0 / val_m0


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

    return param


###### system specific parameter functions
def bulk(bank, material, new_gamma_0=None):
    """Get bulk parameters of a specified material."""
    df_pars = load_params(bank)
    dict_pars = df_pars.loc[material].to_dict()
    dict_pars['gamma_0'] = 1/dict_pars.pop('m_c')
    return renormalize_parameters(dict_pars, new_gamma_0)
