import abc
import copy
import json
import os

import numpy as np
import scipy.linalg as la
# from IPython import display
import sympy

import kwant.continuum

from .misc import spin_matrices, rotate, prettify
from .symbols import momentum
from . import parameters


# Read the cache
def _load_cache():
    """Load cached models.

    File semicon/cache.json should be created on package build.
    """
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    fname = os.path.join(BASE_DIR, 'model_cache.json')
    with open(fname) as f:
        models_cache = json.load(f)
    return models_cache

_models_cache = _load_cache()


def validate_coords(coords):
    """Validate coords in the same way it happens in kwant.continuum."""
    coords = list(coords)
    if coords != sorted(coords):
        raise ValueError("The argument 'coords' must be sorted.")
    if any(c not in 'xyz' for c in coords):
        raise ValueError("The argument 'coords' may only contain "
                         "'x', 'y', or 'z'.")
    return coords


class Model(metaclass=abc.ABCMeta):
    """Simple continuum model.

    Parameters
    ----------
    hamiltonian : sympy.Expr or sympy.Matrix
        Corresponding Hamiltonian.

    Attributes
    ----------
    hamiltonian : str, sympy.Expr or sympy.Matrix
    spin_operators : spin_operators, a 3D tensor
    spins : sequence of spins, alternative to spin_operators
    locals : dict or None, to be passed to kwant.continuum.sympify if
             hamiltonian is string

    Methods
    -------
    rotate : rotate model, see documentation of the method
    prettify : prettify model, see documentation of the meth
    """
    def __init__(self, hamiltonian, spin_operators=None, spins=None,
                 locals=None):
        if isinstance(hamiltonian, str):
            hamiltonian = kwant.continuum.sympify(hamiltonian, locals=locals)
        elif locals is not None:
            raise ValueError('locals can be not None only when hamiltonian is '
                             'of type string.')

        if (spin_operators is not None) and (spins is not None):
            raise ValueError(
                '"spin_operators" and "spins" are mutually exclusive'
            )
        elif spins is not None:
            spin_operators = self.spin_operators(spins)

        if spin_operators is not None:
            expected_shape = (3, *hamiltonian.shape)
            if spin_operators.shape != expected_shape:
                raise ValueError("Shape of spin operators is expected to "
                                 "be {}".format(expected_shape))

        self.hamiltonian = hamiltonian
        self.spin_operators = spin_operators

    def rotate(self, R, act_on=momentum, act_on_spin=True):
        spin_operators = self.spin_operators if act_on_spin else None
        hamiltonian = rotate(self.hamiltonian, R=R, act_on=act_on,
                             spin_operators=spin_operators)

        output = copy.deepcopy(self)
        output.hamiltonian = hamiltonian
        return output

    def prettify(self, decimals=None, zero_atol=None, nsimplify=False):
        hamiltonian = prettify(self.hamiltonian, decimals=decimals,
                               zero_atol=zero_atol, nsimplify=nsimplify)

        output = copy.deepcopy(self)
        output.hamiltonian = hamiltonian
        return output

    @staticmethod
    def spin_operators(spins):
        operators = []
        for s in np.atleast_1d(spins):
            # Explicit if clause seems more clear than oneliner with np.sign
            # spin_matrices: float -> tupple of three spin operators (x, y, z)
            if s > 0:
                operators.append(spin_matrices(s))
            else:
                operators.append(-spin_matrices(-s))

        operators = [
            la.block_diag(*[p[i] for p in operators]) for i in range(3)
        ]

        return np.array(operators)




class BandModel(Model):
    """Basic band-aware model."""

    def __init__(self, bands, components):

        bands = np.atleast_1d(bands)
        components = np.atleast_1d(components)

        # Validate input arguments
        if not all(band in self._allowed_bands for band in bands):
            raise ValueError("Please provide valid bands. Allowed"
                             "bands are {}".format(self._allowed_bands))

        if not all(c in self._allowed_components for c in components):
            raise ValueError("Please provide valid components. Allowed"
                             "components are {}".format(self._allowed_components))

        # If everything is good we proceed with assigning the input arguments
        self.bands = bands
        self.components = components

        # Now we can build hamiltonian and the spin operators
        hamiltonian = self._build_hamiltonian()

        # Finally, we call base constructor to
        Model.__init__(self, hamiltonian=hamiltonian, spins=self.spins)

    @abc.abstractmethod
    def parameters(self, material, databank):
        pass

    @abc.abstractmethod
    def _build_hamiltonian(self):
        pass

    @property
    @abc.abstractmethod
    def _allowed_bands():
        """Sequence of allowed model bands."""
        pass

    @property
    @abc.abstractmethod
    def _allowed_components():
        """Sequence of allowed model components."""
        pass

    @property
    @abc.abstractmethod
    def _band_spins():
        """Mapping: band name -> spin."""
        pass

    @property
    def spins(self):
        spins = [self._band_spins[band] for band in self.bands]
        return spins



class ZincBlende(BandModel):
    """Model for ZincBlende crystals."""

    _allowed_components = ('foreman', 'zeeman')
    _allowed_bands = ('gamma_6c', 'gamma_8v', 'gamma_7v')

    _band_spins = {
        'gamma_6c': 1/2,
        'gamma_8v': 3/2,
        'gamma_7v': 1/2
    }

    _band_indices = {
        'gamma_6c': [0, 1],
        'gamma_8v': [2, 3, 4, 5],
        'gamma_7v': [6, 7]
    }

    _varied_parameters = ['E_0', 'E_v', 'Delta_0', 'P', 'kappa', 'g_c', 'q',
                          'gamma_0', 'gamma_1', 'gamma_2', 'gamma_3']

    def __init__(self, bands=('gamma_6c', 'gamma_8v', 'gamma_7v'),
                 components=('foreman',), parameter_coords=None,
                 default_databank=None):
        self._parameter_coords = parameter_coords

        if isinstance(default_databank, str):
            self.default_databank = parameters.DataBank(default_databank)
        else:
            self.default_databank = default_databank

        BandModel.__init__(self, bands=bands, components=components)

    def _build_hamiltonian(self):
        # return foreman(self._parameter_coords, self.components, self.bands)
        if self._parameter_coords is not None:
            self._parameter_coords = validate_coords(self._parameter_coords)
            str_coords = '({})'.format(", ".join(self._parameter_coords))
            subs = {v: v + str_coords for v in self._varied_parameters}
        else:
            subs = {}

        hamiltonian_components = [
            kwant.continuum.sympify(_models_cache[c], locals=subs)
            for c in self.components
        ]

        hamiltonian = sympy.ImmutableMatrix(sympy.MatAdd(*hamiltonian_components))

        indices = []
        for band in self.bands:
            indices += self._band_indices[band]

        return hamiltonian[:, indices][indices, :]

    def parameters(self, material, databank=None, valence_band_offset=0):
        if databank is None:
            if self.default_databank is not None:
                databank = self.default_databank
            else:
                raise ValueError("No databank provided.")

        output = parameters.ZincBlendeParameters(
            name=material,
            bands=self.bands,
            parameters=databank[material],
            valence_band_offset=valence_band_offset,
        )

        output.update(parameters.constants)
        return output
