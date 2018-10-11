import abc
import copy
import os
import re
import yaml

from collections import UserDict

import numpy as np
from scipy.constants import physical_constants as phys_const

import kwant


# General constants and globals
constants = {
    'm_0': phys_const['electron mass energy equivalent in MeV'][0] * 1e6,
    'phi_0': 2 * phys_const['mag. flux quantum'][0] * (1e9)**2,
    'mu_B': phys_const['Bohr magneton in eV/T'][0],
    'hbar': phys_const['Planck constant over 2 pi times c in MeV fm'][0]
}

taa = constants['hbar']**2 / 2 / constants['m_0']

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATABANK_DIR = os.path.join(BASE_DIR, 'databank')



# Find available databanks
def _find_available_databanks():
    names = []
    for fname in os.listdir(DATABANK_DIR):
        match = re.match("^bank_(.+)\.yml$", fname)
        if match is not None:
            names.append(match.group(1))
    return names


_banks_names = _find_available_databanks()


class DataBank(UserDict):
    """Data bank of effective parameters."""

    def __init__(self, bank):
        UserDict.__init__(self)

        if not isinstance(bank, str):
            raise NotImplementedError("Bank must be a string (for now).")

        if (isinstance(bank, str)) and (bank not in _banks_names):
            msg = "Unkown bankname. Possible options are {}"
            raise ValueError(msg.format(_banks_names))

        self.name = bank
        self.fpath = fpath = os.path.join(DATABANK_DIR, 'bank_' + bank + '.yml')

        with open(fpath, 'r') as f:
            for name, data in yaml.load(f.read()).items():
                self.data[name] = data['parameters']

    def __str__(self):
        output = "Databank:\n"
        output += "    bank name: {}\n".format(self.name)
        output += "    bank path: {}\n".format(self.fpath)
        output += "    materials: " + ", ".join(list(self))
        return output


class Parameters(UserDict):
    """Base parameters class."""
    pass


class BandParameters(Parameters):
    """Basic band-aware parameters class."""

    _loaded_databanks = {}

    def __init__(self, bands, components, material, databank_name):
        Parameters.__init__(self)
        self.bands = bands
        self.components = components
        self.material = material

        if databank_name not in self._loaded_databanks:
            new_databank = DataBank(databank_name)
            self._loaded_databanks[databank_name] = new_databank

        bare = self._calculate_bare(self._loaded_databanks[databank_name])
        self.data.update(bare)

    @property
    @abc.abstractmethod
    def _renormalization_rules(self):
        pass

    def _calculate_bare(self, databank):
        parameters = databank[self.material].copy()
        parameters['gamma_0'] = 1 / parameters.pop('m_c')
        parameters['E_v'] = 0

        renormalizations = self._renormalization_rules

        bare_parameters = parameters.copy()
        for parameter_name in set(parameters) & set(renormalizations):
            rules = renormalizations[parameter_name]
            for band_name in rules:
                # if band not present we can continue
                if band_name not in self.bands:
                    continue
                # otherwise we undo the lowdin transformation
                f = kwant.continuum.lambdify(rules[band_name])
                kwargs = {'T': taa}
                for name in set(f.__code__.co_varnames) - {'T'}:
                    try:
                        kwargs[name] = parameters[name]
                    except KeyError:
                        raise ValueError(
                            "Cannot compute bare value of {}. Parameter "
                            "{} is unkown. Please update databank first."\
                            .format(parameter_name, name)
                        )
                modifier = f(**kwargs)
                bare_parameters[parameter_name] -= modifier

        return {**bare_parameters, **constants}


class ZincBlendeParameters(BandParameters):

    @property
    def _renormalization_rules(self):
        rules = {
            'gamma_0': {
                'gamma_8v': "(2 / 3) * (1 / T) * P**2 / E_0",
                'gamma_7v': "(1 / 3) * (1 / T) * P**2 / (E_0 + Delta_0)"
            },
            'g_c': {
                'gamma_8v': "-(2 / 3) * (1 / T) * P**2 / E_0",
                'gamma_7v': "(2 / 3) * (1 / T) * P**2 / (E_0 + Delta_0)"
            },
            'gamma_1': {
                'gamma_6c': "(1 / 3) * (1 / T) * P**2 / E_0"
            },
            'gamma_2': {
                'gamma_6c': "(1 / 6) * (1 / T) * P**2 / E_0"
            },
            'gamma_3': {
                'gamma_6c': "(1 / 6) * (1 / T) * P**2 / E_0"
            },
            'kappa': {
                'gamma_6c': "(1 / 6) * (1 / T) * P**2 / E_0"
            }
        }
        return rules

    def renormalize(self, new_gamma_0=None, new_P=None):
        if (new_gamma_0 is not None) and (new_P is not None):
            raise ValueError("'new_gamma_0' and 'new_P' are mutually exclusive.")

        output = copy.deepcopy(self)

        f1 = (output['E_0'] * (output['E_0'] + output['Delta_0']))
        f2 = output['E_0'] + (2/3) * output['Delta_0']
        factor = taa * ( f1 / f2 )

        if new_gamma_0 is not None:
            P2 = output['P']**2 + (output['gamma_0'] - new_gamma_0) * factor
            output['P'] = np.sqrt(P2)
            output['gamma_0'] = new_gamma_0
            return output

        if new_P is not None:
            g0 = output['gamma_0'] + (output['P']**2 - new_P**2) / factor
            output['gamma_0'] = g0
            output['P'] = new_P
            return output
