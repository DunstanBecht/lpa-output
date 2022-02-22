#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module analyze.
"""

from lpa.output import analyze
import numpy as np

"""
The following line loads the data of a simulation output. Each piece of
data is stored in a dictionary and can be accessed with a key that can
be found in the documentation of the function output_data().
"""
outdat = analyze.output_data(
    '10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output',
    impdir='data',
)
print(outdat['stm'])
print()

"""
The following line exports in the folder fits/ the fitting data of a
simulation output that can be found in data/. The figures will be
exported in SVG format.
"""
analyze.export(
    '10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output',
    impdir='data',
    expdir='fits',
    figttl=r"10 RDD $ \left( d = 5 \times 10^{-5} \mathrm{nm^{-2}} \right) $",
    fmtout='svg',
    fmtfit='svg',
)

"""
The following line exports in the folder fits/ the fitting data of a
simulation output that can be found in data/. Only the harmonics 1 and
3 will be considered and fitted. The models will be fitted to the
module of the Fourier transform.
"""
analyze.export(
    'rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_PBC1_S0_output.dat',
    impdir='data',
    expdir='fits',
    figttl=r"RDD $ \left( d = 5 \times 10^{-5} \mathrm{nm^{-2}} \right) $",
    frrprt=np.absolute,
    j=np.array((1, 3)),
)

input("OK")
