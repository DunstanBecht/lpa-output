#!/usr/bin/env python
# coding: utf-8

"""
Tools for the analysis of X-ray diffraction simulation output.

After the simulation, the output files containing the cos and sin
coefficients of Fourier analysis can be studied with this package.
The different models are then fitted to the simulation results and for
each fit we can calculate the density and outer cut-off radius
predicted by the model.
"""

__author__ = "Dunstan Becht"
__version__ = "0.9.6"

import os
import sys
import numpy as np
from typing import Union, Optional, NewType, Any

if sys.version_info[0]>=3 and sys.version_info[1]>=9:
    from collections.abc import Callable
else:
    from typing import Callable

if sys.version_info[0]>=3 and sys.version_info[1]>=8:
    from typing import get_args, get_origin
else:
    get_args = lambda x: x.__args__
    get_origin = lambda x: x.__origin__ if hasattr(x, '__origin__') else x

from beartype import beartype

# scalar and vectors
Scalar = Union[int, np.integer, float, np.floating]
Vector = np.ndarray # shape: (n,)

# sets
ScalarList = np.ndarray # shape: (...,)
ScalarListList = Union[list, np.ndarray]
VectorList = np.ndarray # shape: (..., n)

# functions
ModelFunction = Callable[
    [Scalar, Scalar, dict, int, ScalarList],
    ScalarList
]

def getkwa(
    key: str,
    kwa: dict,
    typ: type,
    dft: Any,
) -> Any:
    """
    Retrieve the keyword argument and check its type.

    Input:
        key (str): keyword argument name
        kwa (dict): keyword arguments dictionary
        typ (type): keyword argument type
        dft (Any): keyword argument default value

    Output:
        val (Any): value of the keyword argument
    """
    val = kwa.pop(key, dft) # get the argument or its default value
    if get_origin(typ) == Union: # if the type is a union of types
        typ = get_args(typ) # get the accepted types
    if not isinstance(val, typ): # if the type is not correct
        msg = (f"keyword argument {key!r} must be of type {typ.__name__!r} "
               f"but {val!r} is of type {type(val).__name__!r}")
        raise TypeError(msg)
    return val

def endkwa(
    kwa: dict,
) -> None:
    """
    Check that no bad keyword arguments have been passed.

    Input:
        kwa (dict): keyword arguments dictionary
    """
    if len(kwa) > 0: # if there are still keyword arguments after recovery
        msg = "unexpected keyword argument: {!r}={!r}".format(*kwa.popitem())
        raise TypeError(msg)
