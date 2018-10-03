import abc
import copy

from collections import UserDict

import numpy as np

from .misc import spin_matrices, rotate, prettify
from .symbols import momentum

class Parameters(UserDict):
    """Base parameters class

    Note: inheriting from "UserDict" makes "Parameters" an abstract class.
    """
    pass


class Model(metaclass=abc.ABCMeta):
    """Base class for all models."""
    pass


class BulkModel(Model):
    """Simple bulk model.

    Parameters
    ----------
    hamiltonian : sympy.Expr or sympy.Matrix
        Corresponding Hamiltonian.
    parameters : dict-like (optional)
        Hamiltonian parameters.

    Attributes
    ----------
    hamiltonian : sympy.Expr or sympy.Matrix
    spins : sequence of integers or half-integers

    Properties
    ----------
    parameters : Parameters instance

    Methods
    -------
    rotate : rotate model, see documentation of the method
    """
    def __init__(self, hamiltonian, parameters=None, spins=None):
        self.hamiltonian = hamiltonian
        self.parameters = parameters if parameters is not None else {}
        self.spins = spins = np.atleast_1d(spins)

        if spins is not None:
            d = sum([2 * s + 1 for s in spins])
            if int(d) != hamiltonian.shape[0]:
                raise ValueError("Dimension of spins space doesn't match "
                                 "dimension of the Hamiltonian.")

    @property
    def spin_operators(self):
        """Return spin operators of defined spins."""

        if self.spins is None:
            return None

        operators = []
        for s in self.spins:
            # Explicit if clause seems more clear than oneliner with np.sign
            if s > 0:
                operators.append(spin_matrices(s))
            else:
                operators.append(-spin_matrices(-s))

        operators = [np.bmat([p[i] for p in operators]) for i in range(3)]

        return np.array(operators)

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


# class BandModel(BulkModel):

#     def __init__(self, bands, components):

#         # I am not sure if this is not too specific and should be handled
#         # by a subclass constructor
#         self.bands = bands
#         self.components = components
#         self.coords = coords

#         # self._hamiltonian = ...


#     @property
#     @abc.abstractmethod
#     def hamiltonian(self):
#         pass

#     # @hamiltonian.setter
#     # def hamiltonian(self, val):
#     #     self._hamiltonian = val

#     @property
#     @abc.abstractmethod
#     def spin_operators(self):
#         pass

#     @abc.abstractmethod
#     def parameters(self, parameter):
#         pass

#     # def rotate(self, R):
#     #     # I want this method to be the same for all subclasses
#     #     # however I am not sure how shall I handle modification of Hamiltonians
#     #     self.hamiltonian = rotate(R, self.hamiltonian, self.spin_operators)




# class Zincblende(Model):
#     """Model for Zincblende crystals."""

#     def __init__(self, bands, components, coords=None):

#         # Do some prechecks on parameters

#         # Call Super constructor
#         Model.__init__(self, bands, components, coords)
#         pass

#     @property
#     def hamiltonian(self):
#         """Return model Hamiltonian"""

#         # Code that generate the actual Hamiltonian
#         pass
