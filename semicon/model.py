import abc
from collections import UserDict

from .tools import rotate


class Parameters(UserDict):
    """Base parameters class

    Note: inheriting from "UserDict" makes "Parameters" an abstract class.
    """
    pass



class Model(metaclass=abc.ABCMeta):
    """Base class for all models.

    Attributes
    ----------
    hamiltonian : sympy.Expr or sympy.Matrix
        Hamiltonian corresponding to the specified model.

    spin_operators : sequence of matrices
        Corresponding spin operators.

    Methods
    -------
    parameters : Parameters instance (or its subclass)
        Corresponding parameters

    rotate : Model
        Apply rotation to the model
    """
    def __init__(self, bands, components, coords=None):

        # I am not sure if this is not too specific and should be handled
        # by a subclass constructor
        self.bands = bands
        self.components = components
        self.coords = coords

        # self._hamiltonian = ...


    @property
    @abc.abstractmethod
    def hamiltonian(self):
        pass

    # @hamiltonian.setter
    # def hamiltonian(self, val):
    #     self._hamiltonian = val

    @property
    @abc.abstractmethod
    def spin_operators(self):
        pass

    @abc.abstractmethod
    def parameters(self, parameter):
        pass

    # def rotate(self, R):
    #     # I want this method to be the same for all subclasses
    #     # however I am not sure how shall I handle modification of Hamiltonians
    #     self.hamiltonian = rotate(R, self.hamiltonian, self.spin_operators)




class Zincblende(Model):
    """Model for Zincblende crystals."""

    def __init__(self, bands, components, coords=None):

        # Do some prechecks on parameters
        ...

        # Call Super constructor
        Model.__init__(self, bands, components, coords)

    @property
    def hamiltonian(self):
        """Return model Hamiltonian"""

        # Code that generate the actual Hamiltonian
        ...

        return ham
