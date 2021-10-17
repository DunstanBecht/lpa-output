#!/usr/bin/env python
# coding: utf-8

"""
Available models.
"""

import scipy
import scipy.integrate
from . import *

@beartype
def GUW1(
    p: ScalarList,
    o: dict,
    j: int,
    l: ScalarList,
) -> ScalarList:
    """
    Return the Fourier amplitudes calculated with the model GUW1.

    The model is described in: I. Groma, T. Ungár, and M. Wilkens.
    “Asymmetric X-ray line broadening of plastically deformed crystals.
    I. Theory”. In: Journal of Applied Crystallography (1988). ISSN:
    0021-8898.

    Input:
        p (ScalarList): optimization parameters, contains:
            d (Scalar): density of dislocations [nm^-2]
            r (Scalar): outer cut-off radius [nm]
            f (Scalar): fluctuation oof the density [1]
            R (Scalar): R0 [nm]
        o (dict): output data dictionary
        j (int): selected harmonic
        l (ScalarList): Fourier variable [nm]

    Output:
        a (ScalarList): Fourier amplitudes
    """
    d, r, f, R = p
    i = o['index'][j] # index of the harmonic in c
    jgb, jg2, b2, C = o['jgb'][i], o['jg2'][i], o['b2'], o['C']
    k = np.pi/2*jg2*b2*C*d
    D = k * (np.log(l/r)+f/2*k*l**2*np.log(l/R)*2)
    return np.exp(l**2*D)

@beartype
def GUW2(
    p: ScalarList,
    o: dict,
    j: int,
    l: ScalarList,
) -> ScalarList:
    """
    Return the Fourier amplitudes calculated with the model GUW2.

    The model is described in: I. Groma, T. Ungár, and M. Wilkens.
    “Asymmetric X-ray line broadening of plastically deformed crystals.
    I. Theory”. In: Journal of Applied Crystallography (1988). ISSN:
    0021-8898.

    Input:
        p (ScalarList): optimization parameters, contains:
            d (Scalar): density of dislocations [nm^-2]
            r (Scalar): outer cut-off radius [nm]
        o (dict): output data dictionary
        j (int): selected harmonic
        l (ScalarList): Fourier variable [nm]

    Output:
        a (ScalarList): Fourier amplitudes
    """
    d, r = p
    i = o['index'][j] # index of the harmonic in c
    jgb, jg2, b2, C = o['jgb'][i], o['jg2'][i], o['b2'], o['C']
    k = np.pi/2*jg2*b2*C*d
    D = k * (np.log(l)-np.log(jgb*r))
    return np.exp(l**2*D)

@beartype
def f(
    e: Scalar,
) -> Scalar:
    """
    Intermediary function in the calculation of the model W1.

    Input:
        e (Scalar): Fourier variable divided by the outer cut-off radius [1]

    Output:
        r (Scalar): value of f on e
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
        r = 512/90/np.pi/e - (11/24+np.log(2)*e/4)/e2
    return r

vf = np.vectorize(f)

@beartype
def W1(
    p: ScalarList,
    o: dict,
    j: int,
    l: ScalarList,
) -> ScalarList:
    """
    Return the Fourier amplitudes calculated with the model W1.

    The model is described in: M. Wilkens. Fundamental aspects of
    dislocation theory. Ed. by J. A. Simmons, R. de Wit, and R.
    Bullough. Vol. 2. U.S. National Bureau of Standards, 1970, pp.
    1195–1221.

    Input:
        p (ScalarList): optimization parameters contains:
            d (Scalar): density of dislocations [nm^-2]
            r (Scalar): outer cut-off radius [nm]
        o (dict): output data dictionary
        j (int): selected harmonic
        l (ScalarList): Fourier variable [nm]

    Output:
        a (ScalarList): Fourier amplitudes
    """
    d, r = p
    i = o['index'][j] # index of the harmonic in c
    jg2, b2, C = o['jg2'][i], o['b2'], o['C']
    k = np.pi/2*jg2*b2*C*d
    D = -k * vf(np.exp(1/4)*l/2/r)
    AD = np.exp(l**2*D)
    AS = 1
    return AD * AS

@beartype
def W2(
    p: ScalarList,
    o: dict,
    j: int,
    l: ScalarList,
) -> ScalarList:
    """
    Return the Fourier amplitudes calculated with the model W2.

    The model is described in: J.-D. Kamminga and R. Delhez.
    “Calculation of diffraction line profiles for structures with
    dislocations”. In: Materials Science Forum (2001). ISSN: 1662-9752.

    Input:
        p (ScalarList): optimization parameters contains:
            d (Scalar): density of dislocations [nm^-2]
            r (Scalar): outer cut-off radius [nm]
        o (dict): output data dictionary
        j (int): selected harmonic
        l (ScalarList): Fourier variable [nm]

    Output:
        a (ScalarList): Fourier amplitudes
    """
    d, r = p
    i = o['index'][j] # index of the harmonic in c
    jgb, jg2, b2, C = o['jgb'][i], o['jg2'][i], o['b2'], o['C']
    z, g = o['z'], o['g']
    nz, ng = np.linalg.norm(z), np.linalg.norm(g)
    s = np.sqrt(1-(np.dot(z, g)/(nz*ng))**2)
    k = np.pi/2*jg2*b2*C*d
    D = k * (np.log(l)-np.log(jgb*r)-2*np.log(2)+1/3+np.log(s*np.abs(jgb)))
    return np.exp(l**2*D)
