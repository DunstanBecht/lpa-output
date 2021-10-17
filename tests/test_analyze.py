#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module analyze.
"""

from lpa.output import analyze
import numpy as np

# test output_data
outdat = analyze.output_data(
    '10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output',
    impdir='data',
)
print(outdat['stm'])
print()

# test export
analyze.export(
    '10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output',
    impdir='data',
    expdir='fits',
    figttl=r"10 RDD $ \left( d = 5 \times 10^{-5} \mathrm{nm^{-2}} \right) $",
    fmtout='svg',
    fmtfit='svg',
)
analyze.export(
    'rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_PBC1_S0_output.dat',
    impdir='data',
    expdir='fits',
    figttl=r"RDD $ \left( d = 5 \times 10^{-5} \mathrm{nm^{-2}} \right) $",
    frrprt=np.absolute,
    j=np.array((1, 3)),
)


input("OK")
