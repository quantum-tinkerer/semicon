# Function defined in this file serve miscellaneous purposes.
# See comments for each group of functions for more details.


import warnings

import numpy as np
import scipy.linalg as la
from collections import defaultdict

from scipy.interpolate import interp1d

try:
    from scipy.spatial.transform import Rotation
    rotation_functionality_available = True
except ImportError:
    warnings.warn("Rotation functionality not availble unless "
                  "scipy 1.2 or greater is installed", RuntimeWarning)
    rotation_functionality_available = False

import sympy
import kwant.continuum

from .symbols import momentum


def spin_matrices(s):
    """Construct spin-s matrices for any spin."""
    d = np.round(2*s + 1)
    if not np.isclose(d, int(d)):
        raise ValueError("Argument 's' must be integer or half integer.")
    d = int(d)
    Sz = 1/2 * np.diag(np.arange(d - 1, -d, -2))
    # first diagonal for general s from en.wikipedia.org/wiki/Spin_(physics)
    diag = [(1 / 2) * np.sqrt(2 * i * (s + 1) - i * (i + 1)) for i in np.arange(1, d)]
    Sx = np.diag(diag, k=1) + np.diag(diag, k=-1)
    Sy = -1j * np.diag(diag, k=1) + 1j * np.diag(diag, k=-1)
    return np.array([Sx, Sy, Sz])


def _prettify_term(expr, decimals=None, zero_atol=None, nsimplify=False):
    terms = [(k, v) for k, v in monomials(expr).items()]

    output = []
    for k, v in terms:
        v = complex(v)

        # numerical rounding to given precision
        if decimals is not None:
            v = np.round(v, decimals)

        # check if v is zero up to desired atol
        if zero_atol is not None:
            v_real = 0 if np.isclose(v.real, 0, atol=zero_atol) else v.real
            v_imag = 0 if np.isclose(v.imag, 0, atol=zero_atol) else v.imag
            v = v_real + 1j * v_imag

        # sympy nsimplify (subs sqrt(3) and similar in place of floats)
        if nsimplify:
            v = sympy.nsimplify(sympy.sympify(v))

        output.append(k * v)
    return sympy.Add(*output)


def prettify(expression, decimals=None, zero_atol=None, nsimplify=False):
    """Prettify SymPy expression.

    1. cast to monomials: symbols -> numerical factor
    2. if decimals is not None: use np.round on numerical factor
    3. if zero_atol is not None: check np.isclose(x, 0, atol=zero_atol)
       to check if number is zero
    4. if nsimplify is True: use sympy.nsimplify
    """
    if not isinstance(expression, sympy.matrices.MatrixBase):
        return _prettify_term(expression, decimals, zero_atol, nsimplify)

    output = sympy.zeros(*expression.shape)
    for (i, j), expr in np.ndenumerate(expression):
        output[i, j] = _prettify_term(expr, decimals, zero_atol, nsimplify)

    return output


def sympy_to_numpy(arr, dtype=complex):
    if isinstance(arr, np.ndarray):
        return arr
    else:
        return np.array(arr.tolist(), dtype=dtype)


def _validate_rotation_matrix(R):
    if not R.shape[0] == R.shape[1] == 3:
        raise ValueError("Rotation matrix is expected to be 3x3 unitary Matrix")
    if isinstance(R, np.ndarray):
        det = la.det(R)
    elif isinstance(R, sympy.matrices.MatrixBase):
        det = R.det()
    else:
        raise ValueError("rotation matrix should be defined as np.array or "
                         "sympy.Matrix")

    if not np.allclose(float(det), 1):
        raise ValueError("Determinant of rotation matrix must be 1.")


def _basis_rotation(R, spin_operators):
    R = sympy_to_numpy(R, dtype=float)
    n = Rotation.from_dcm(R).as_rotvec()
    spin_operators = [sympy_to_numpy(s) for s in spin_operators]
    ns = np.sum([ni * si for (ni, si) in zip(n, spin_operators)], axis=0)
    return la.expm(1j * ns)


def rotate(expr, R, act_on=momentum, spin_operators=None):
    if not rotation_functionality_available:
        raise RuntimeError("Rotation functionality not availble. Please, "
                           "install scipy 1.2 or greater.")

    _validate_rotation_matrix(R)
    rotation_subs = lambda R, v: {cprime: c for (cprime, c) in zip(v, R @ v)}

    subs = {}
    for row in np.atleast_2d(act_on):
        new_subs = rotation_subs(R, sympy.Matrix(row))
        subs.update(new_subs)

    # "simultaneous" flag is very important here
    # note that SymPy takes it as **kwargs and there
    # is no validation for typos!!!
    expr = expr.subs(subs, simultaneous=True)

    if spin_operators is not None:
        U = _basis_rotation(R, spin_operators)
        expr = U @ expr @ U.transpose().conjugate()

    return expr.expand()



# Function defined in this section come from "kwant.continuum" module
# of Kwant and are currently a part of a non-public API.
# To avoid breakage with future releases, they are defined here.

def make_commutative(expr, *symbols):
    """Make sure that specified symbols are defined as commutative.

    Parameters
    ----------
    expr: sympy.Expr or sympy.Matrix
    symbols: sequace of symbols
        Set of symbols that are requiered to be commutative. It doesn't matter
        of symbol is provided as commutative or not.

    Returns
    -------
    input expression with all specified symbols changed to commutative.
    """
    names = [s.name if not isinstance(s, str) else s for s in symbols]
    symbols = [sympy.Symbol(name, commutative=False) for name in names]
    expr = expr.subs({s: sympy.Symbol(s.name) for s in symbols})
    return expr


def monomials(expr, gens=None):
    """Parse ``expr`` into monomials in the symbols in ``gens``.

    Parameters
    ----------
    expr: sympy.Expr or sympy.Matrix
        Sympy expression to be parsed into monomials.
    gens: sequence of sympy.Symbol objects or strings (optional)
        Generators of monomials. If unset it will default to all
        symbols used in ``expr``.

    Returns
    -------
    dictionary (generator: monomial)

    Example
    -------
        >>> expr = kwant.continuum.sympify("A * (x**2 + y) + B * x + C")
        >>> monomials(expr, gens=('x', 'y'))
        {1: C, x: B, x**2: A, y: A}
    """
    if gens is None:
        gens = expr.atoms(sympy.Symbol)
    else:
        gens = [kwant.continuum.sympify(g) for g in gens]

    if not isinstance(expr, sympy.MatrixBase):
        return _expression_monomials(expr, gens)
    else:
        output = defaultdict(lambda: sympy.zeros(*expr.shape))
        for (i, j), e in np.ndenumerate(expr):
            mons = _expression_monomials(e, gens)
            for key, val in mons.items():
                output[key][i, j] += val
        return dict(output)


def _expression_monomials(expr, gens):
    """Parse ``expr`` into monomials in the symbols in ``gens``.

    Parameters
    ----------
    expr: sympy.Expr
        Sympy expr to be parsed.
    gens: sequence of sympy.Symbol
        Generators of monomials.

    Returns
    -------
    dictionary (generator: monomial)
    """
    expr = sympy.expand(expr)
    output = defaultdict(lambda: sympy.Integer(0))
    for summand in expr.as_ordered_terms():
        key = []
        val = []
        for factor in summand.as_ordered_factors():
            symbol, exponent = factor.as_base_exp()
            if symbol in gens:
                key.append(factor)
            else:
                val.append(factor)
        output[sympy.Mul(*key)] += sympy.Mul(*val)

    return dict(output)



### Helper functions, to be replaced with something better...
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
