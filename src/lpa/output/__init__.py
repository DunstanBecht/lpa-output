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
__version__ = "0.8.5"

import os
import sys
import numpy as np
from typing import Union, Optional, NewType, Any
from beartype import beartype

if sys.version_info[0]>=3 and sys.version_info[1]>=9:
    from collections.abc import Callable
    List = list
    Tuple = tuple
else:
    from typing import Callable
    from typing import List
    from typing import Tuple

# scalar and vectors
Scalar = Union[int, float, np.intc]
Vector = np.ndarray # shape: (n,)
# sets
ScalarList = np.ndarray # shape: (...,)
ScalarListList = Union[list, np.ndarray]
VectorList = np.ndarray # shape: (..., n)
# functions
ModelFunction = Callable[
    [Scalar, Scalar, dict, int, ScalarList],
    [ScalarList]
]
