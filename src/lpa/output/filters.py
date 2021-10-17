#!/usr/bin/env python
# coding: utf-8

"""
Filters to delimit the range of the Fourier variable.
"""

from . import *

@beartype
def f0(
    a: ScalarList,
) -> int:
    """
    Return the index of the first negative value.

    Input:
        a: Fourier amplitudes for a given harmonic

    Output:
        i0: index of the first negative value
    """
    i0 = 0
    while i0<len(a) and a[i0]>0:
        i0 += 1
    return i0

@beartype
def f1(
    a: ScalarList,
) -> int:
    """
    Return the index from which the noise begins.

    Input:
        a: Fourier amplitudes for a given harmonic

    Output:
        i1: index at which the noise starts
    """
    i0 = f0(a)
    i1 = 1
    while i1<i0 and a[i1-1]>a[i1]:
        i1 += 1
    return i1

@beartype
def f2(
    a: ScalarList,
    l: ScalarList,
) -> int:
    """
    Return the index from which ln(A(L))/L^2 is no longer linear.

    Input:
        l: Fourier variable
        a: Fourier amplitudes for a given harmonic

    Output:
        i2: index that marks the end of the linear part
    """
    i0 = f0(a)
    a, l = a[:i0], l[:i0]
    # convert
    y = np.log(a)/l**2
    x = np.log(l)
    # filter
    return f2_xy(y, x)

@beartype
def f2_xy(
    y: ScalarList,
    x: ScalarList,
) -> int:
    """
    Return the index from which y is no longer linear.

    Input:
        y: y variable
        x: x variable

    Output:
        i2: index that marks the end of the linear part
    """
    w = np.ones(len(x))
    xw = np.stack((x, w)).T
    if len(x)<3:
        raise ValueError("not enough points")
    err = 0
    i = 3
    while i<len(x) and np.abs(err/np.ptp(y[:i]))<0.01:
        i += 1
        err = np.sqrt(np.linalg.lstsq(xw[:i], y[:i], rcond=None)[1][0]/(i-2))
    return i
