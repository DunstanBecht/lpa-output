#!/usr/bin/env python
# coding: utf-8

"""
Filters to delimit the range of the Fourier variable.
"""

from . import *

@beartype
def F0(
    a: ScalarList,
) -> int:
    """
    Return the index of the first negative value.

    Input:
        a (ScalarList): Fourier amplitudes for a given harmonic

    Output:
        i0 (int): index of the first negative value
    """
    i0 = 0
    while i0<len(a) and a[i0]>0:
        i0 += 1
    return i0

@beartype
def F1(
    a: ScalarList,
) -> int:
    """
    Return the index from which the noise begins.

    Input:
        a (ScalarList): Fourier amplitudes for a given harmonic

    Output:
        i1 (int): index at which the noise starts
    """
    i0 = F0(a)
    i1 = 1
    while i1<i0 and a[i1-1]>a[i1]:
        i1 += 1
    return min(i1, i0)

@beartype
def F2(
    a: ScalarList,
    l: ScalarList,
) -> int:
    """
    Return the index from which ln(A(L))/L^2 is no longer linear.

    Input:
        l (ScalarList): Fourier variable
        a (ScalarList): Fourier amplitudes for a given harmonic

    Output:
        i2 (int): index that marks the end of the linear part
    """
    i0 = F0(a)
    a, l = a[:i0], l[:i0]
    # convert
    y = np.log(a)/l**2
    x = np.log(l)
    # filter
    i2 = F2_xy(y, x)
    return min(i2, i0)

@beartype
def F2_xy(
    y: ScalarList,
    x: ScalarList,
) -> int:
    """
    Return the index from which y is no longer linear.

    Input:
        y (ScalarList): y variable
        x (ScalarList): x variable

    Output:
        i2 (int): index that marks the end of the linear part
    """
    w = np.ones(len(x))
    xw = np.stack((x, w)).T
    if len(x) < 3:
        return 2
    std = lambda a, b: np.sqrt(np.linalg.lstsq(a, b, rcond=None)[1][0]/(i-1))
    i = 3
    while i<len(x) and std(xw[:i+1], y[:i+1])/np.ptp(y[:i+1])<0.01:
        i += 1
    return i
