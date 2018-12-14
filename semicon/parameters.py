import abc
import os
import re
import yaml

import inspect
from collections import UserDict

import numpy as np
import pandas as pd
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
        match = re.match(r"^bank_(.+)\.yml$", fname)
        if match is not None:
            names.append(match.group(1))
    return names


_banks_names = _find_available_databanks()


class DataBank(UserDict):
    """Data bank of effective parameters."""
    def __init__(self, name):
        # If "name" is one of predefined databank then load it, otherwise
        # check if it is absolute path to existin datafile.
        if name in _banks_names:
            fpath = os.path.join(DATABANK_DIR, 'bank_' + name + '.yml')
        elif os.path.isabs(name):
            fpath = name
        else:
            msg = ("Wrong file name. Please provide a valid path to a file "
                   "or choose from the following predefined banks: {}.")
            raise ValueError(msg.format(_banks_names))

        self.name = name

        # Initialize base class
        UserDict.__init__(self)

        # Read parameters from file
        with open(fpath, 'r') as f:
            for name, data in yaml.load(f.read()).items():
                self.data[name] = data['parameters']

    def __str__(self):
        output = "Databank:\n"
        output += "    bank name: {}\n".format(self.name)
        output += "    materials: " + ", ".join(list(self))
        return output

    def to_dataframe(self):
        return pd.DataFrame(self.data).T


class BareParameters(UserDict):
    """Basic band-aware bare parameters class."""

    def __init__(self, name, bands, parameters, already_bare=False):
        self.name = name
        self.bands = bands

        if not already_bare:
            parameters = self._calculate_bare(parameters.copy())
        else:
            parameters = parameters.copy()

        UserDict.__init__(self, **parameters)

    @property
    @abc.abstractmethod
    def _renormalization_rules(self):
        pass

    def to_effective(self):
        return self._calculate_bare(self.data, reverse=True)

    def _calculate_bare(self, parameters, reverse=False):
        renormalizations = self._renormalization_rules

        bare_parameters = parameters.copy()
        for parameter_name in set(parameters) & set(renormalizations):
            # First we go over all renormalization rules for each parameter
            # and if rule-corresponding band is present in bands we apply it
            rules = renormalizations[parameter_name]
            for band_name in rules:
                # if band not present we can continue
                if band_name not in self.bands:
                    continue
                # otherwise we undo the lowdin transformation
                f = kwant.continuum.lambdify(rules[band_name])
                kwargs = {'T': taa}
                for name in set(inspect.signature(f).parameters) - {'T'}:
                    try:
                        kwargs[name] = parameters[name]
                    except KeyError:
                        raise ValueError(
                            "Cannot compute bare value of {}. Parameter "
                            "{} is unkown. Please update databank first."\
                            .format(parameter_name, name)
                        )
                modifier = f(**kwargs)
                if not reverse:
                    bare_parameters[parameter_name] -= modifier
                else:
                    bare_parameters[parameter_name] += modifier

        return bare_parameters


class ZincBlendeParameters(BareParameters):
    """Parameter class for ZincBlende materials."""

    _renormalization_rules = {
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

    def __init__(self, name, bands, parameters, valence_band_offset=0,
                 already_bare=False):
        parameters = parameters.copy()

        if 'm_c' in parameters:
            parameters['gamma_0'] = 1 / parameters.pop('m_c')

        if 'E_v' in parameters:
            parameters['E_v'] += valence_band_offset
        else:
            parameters['E_v'] = valence_band_offset

        BareParameters.__init__(
            self, name=name, bands=bands, parameters=parameters,
            already_bare=already_bare
        )

    def renormalize(self, new_gamma_0=None, new_P=None):
        if (new_gamma_0 is not None) and (new_P is not None):
            msg = "'new_gamma_0' and 'new_P' are mutually exclusive."
            raise ValueError(msg)

        if 'gamma_6c' not in self.bands:
            msg = "Cannot apply workaround without the electron band."
            raise ValueError(msg)

        if ('gamma_7v' not in self.bands) or ('gamma_8v' not in self.bands):
            msg = "Cannot apply workaround without at least one hole band."
            raise ValueError(msg)

        effective = self.to_effective()

        if new_gamma_0 is not None:
            # First, calculate scaling factor
            factor = 0

            if 'gamma_7v' in self.bands:
                factor += (2 / 3) / effective['E_0']

            if 'gamma_8v' in self.bands:
                factor += (1 / 3) / (effective['E_0'] + effective['Delta_0'])

            # Second, calculate required P
            P2 = (effective['gamma_0'] - new_gamma_0) * (taa / factor)
            new_P = np.sqrt(P2)

        effective['P'] = new_P
        output = ZincBlendeParameters(name=self.name, bands=self.bands,
                                      parameters=effective)

        return output
