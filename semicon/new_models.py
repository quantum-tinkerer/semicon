import abc
import copy

import numpy as np
import scipy.linalg as la
# from IPython import display

from .misc import spin_matrices, rotate, prettify
from .symbols import momentum
from . import new_parameters as parameters



class Model(metaclass=abc.ABCMeta):
    """Simple continuum model.

    Parameters
    ----------
    hamiltonian : sympy.Expr or sympy.Matrix
        Corresponding Hamiltonian.

    Attributes
    ----------
    hamiltonian : sympy.Expr or sympy.Matrix
    spin_operators : spin_operators, a 3D tensor
    spins : sequence of spins, alternative to spin_operators

    Methods
    -------
    rotate : rotate model, see documentation of the method
    prettify : prettify model, see documentation of the meth
    """
    def __init__(self, hamiltonian, spin_operators=None, spins=None):
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

        # Validate input arguments
        if not all(band in self._allowed_bands for band in bands):
            raise ValueError("Please provide valid bands. Allowed"
                             "bands are {}".format(self._allowed_bands))

        if not all(c in self._allowed_components for c in components):
            raise ValueError("Please provide valid components. Allowed"
                             "components are {}".format(self._allowed_bands))

        # If everything is good we proceed with assigning the input arguments
        self.bands = bands
        self.components = components

        # Now we can build hamiltonian and the spin operators
        hamiltonian = self._build_hamiltonian()

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
        """Mapping: allowed band -> band's spin."""
        pass

    @property
    @abc.abstractmethod
    def _allowed_components():
        pass

    @property
    def spins(self):
        spins = [self._allowed_bands[band] for band in self.bands]
        return spins




from .models import foreman

class Zincblende(BandModel):
    """Model for Zincblende crystals."""

    _allowed_components = ('foreman', 'zeeman')
    _allowed_bands = {'gamma_6c': 1/2, 'gamma_8v': 3/2, 'gamma_7v': 1/2}

    def __init__(self, bands, components, parameter_dependence=None,
                 default_databank='lawaetz'):
        self._parameter_dependence = parameter_dependence
        self.default_databank = default_databank

        BandModel.__init__(self, bands=bands, components=components)

    def _build_hamiltonian(self):
        return foreman(self._parameter_dependence, self.components, self.bands)

    def parameters(self, material, databank_name=None):
        if databank_name is None:
            databank_name = self.default_databank

        return parameters.ZincBlendeParameters(
            bands=self.bands,
            components=self.components,
            material=material,
            databank_name=databank_name
        )
