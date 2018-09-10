# Function defined in this file serve miscellaneous purposes.
# See comments for each group of functions for more details.


import numpy as np
import scipy.linalg as la
import xarray as xr
from itertools import product
from collections import Mapping, defaultdict

from scipy.spatial.transform import Rotation

import sympy
import kwant

from .kp_models import symbols


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


def rotate_symbols(expr, R):
    rotation_subs = lambda R, v: {cprime: c for (cprime, c) in zip(v, R @ v)}

    subs1 = rotation_subs(R, sympy.Matrix(symbols.momentum_symbols))
    subs2 = rotation_subs(R, sympy.Matrix(symbols.position_symbols))
    subs3 = rotation_subs(R, sympy.Matrix(symbols.magnetic_symbols))
    subs = {**subs1, **subs2, **subs3}
    # "simultaneous" flag is very important here
    # note that SymPy takes it as **kwargs and there
    # is no validation for typos!!!
    return expr.subs(subs, simultaneous=True).expand()


def basis_rotation(R, spin_operators):
    R = sympy_to_numpy(R, dtype=float)
    n = Rotation.from_dcm(R).as_rotvec()
    spin_operators = [sympy_to_numpy(s) for s in spin_operators]
    ns = np.sum([ni * si for (ni, si) in zip(n, spin_operators)], axis=0)
    return la.expm(1j * ns)



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
