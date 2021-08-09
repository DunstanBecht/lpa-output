#!/usr/bin/env python
# coding: utf-8

"""
Studied models.
"""

import scipy
import scipy.integrate
from . import *

@beartype
def Groma(
    d: Scalar,
    r: Scalar,
    c: dict,
    j: int,
    l: ScalarList,
) -> ScalarList:
    """
    Return the Fourier amplitudes calculated with the Groma model.

    Input:
        d: density of dislocations [nm^-2]
        r: outer cut-off radius [nm]
        c: common quantities dictionary
        j: selected harmonic
        l: Fourier variable [nm]

    Output:
        a: Fourier amplitudes
    """
    i = c['index'][j]
    jgb, jg2, b2, C = c['jgb'][i], c['jg2'][i], c['b2'], c['C']
    k = np.pi/2*jg2*b2*C*d
    D = k * (np.log(l)-np.log(jgb*r))
    return np.exp(l**2*D)

@beartype
def Kamminga(
    d: Scalar,
    r: Scalar,
    c: dict,
    j: int,
    l: ScalarList,
) -> ScalarList:
    """
    Return the Fourier amplitudes calculated with the Kamminga model.

    Input:
        d: density of dislocations [nm^-2]
        r: outer cut-off radius [nm]
        c: common quantities dictionary
        j: selected harmonic
        l: Fourier variable [nm]

    Output:
        a: Fourier amplitudes
    """
    i = c['index'][j]
    jgb, jg2, b2, C = c['jgb'][i], c['jg2'][i], c['b2'], c['C']
    z, g = c['z'], c['g']
    nz, ng = np.linalg.norm(z), np.linalg.norm(g)
    s = np.sqrt(1-(np.dot(z, g)/(nz*ng))**2)
    k = np.pi/2*jg2*b2*C*d
    D = k * (np.log(l)-np.log(jgb*r)-2*np.log(2)+1/3+np.log(s*np.abs(jgb)))
    return np.exp(l**2*D)

@beartype
def f(
    e: Scalar,
) -> Scalar:
    """
    Intermediary function in the calculation of the Wilkens model.

    Input:
        e: value of Fourier variable divided by the outer cut-off radius [1]

    Output:
        r: value of f on eta
    """
    e2 = e**2
    def h(v):
        if v==0:
            return 1
        return np.arcsin(v)/v
    if e < 1:
        r = (7/4 - np.log(e) - np.log(2) + 512/90/np.pi/e
            + 2/np.pi*(1-1/4/e2)*scipy.integrate.quad(h, 0, e)[0]
            - 1/np.pi*(769/180/e+41*e/90+2*e2*e/90)*(1-e2)**(1/2)
            - 1/np.pi*(11/12/e2+7/2+e2/3)*np.arcsin(e) + e2/6)
    else:
        r = 512/90/pi/eta - (11/24+np.log(2)*eta/4)/eta2
    return r

vf = np.vectorize(f)

@beartype
def Wilkens(
    d: Scalar,
    r: Scalar,
    c: dict,
    j: int,
    l: ScalarList,
) -> ScalarList:
    """
    Return the Fourier amplitudes calculated with the Kamminga model.

    Input:
        d: density of dislocations [nm^-2]
        r: outer cut-off radius [nm]
        c: common quantities dictionary
        j: selected harmonic
        l: Fourier variable [nm]

    Output:
        a: Fourier amplitudes
    """
    i = c['index'][j]
    jg2, b2, C = c['jg2'][i], c['b2'], c['C']
    k = np.pi/2*jg2*b2*C*d
    D = -k * vf(np.exp(1/4)*l/2/r)
    AD = np.exp(l**2*D)
    AS = 1
    return AD * AS
