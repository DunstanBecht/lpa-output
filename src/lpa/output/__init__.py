#!/usr/bin/env python
# coding: utf-8

"""
Tools for the analysis of X-ray diffraction simulation output.

After the simulation, the output file containing the cos and sin
coefficients from Fourier analysis can be studied with this package.
The different models are then fitted to the simulation results and for
each fit we can calculate the density predicted by the model.
"""

__author__ = "Dunstan Becht"
__version__ = "0.8.3"

import os
import sys
import numpy as np
from typing import Union, Optional, NewType, Any
try:
    from beartype import beartype
except:
    def beartype(function):
        return function

if sys.version_info[0]>=3 and sys.version_info[1]>=9:
    from collections.abc import Callable
    List = list
    Tuple = tuple
else:
    from typing import Callable
    from typing import List
    from typing import Tuple

Scalar = Union[int, float, np.intc]
Vector = np.ndarray # shape: (n,)
ScalarList = np.ndarray # shape: (...,)
ScalarListList = list
VectorList = np.ndarray # shape: (..., n)
