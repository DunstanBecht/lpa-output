#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module analyze.
"""

from lpa.output import analyze
import numpy as np

stm = '10_rho5e14m-2_square_2000nm_RDD_d5e-4nm-2_edge_PBCR2_S0'
ttl = r"10 RDD $ \left( d = 5 \times 10^{-4} \mathrm{nm^{-2}} \right) $"

# test output_data
o = analyze.output_data(stm, 'data')
print(o['stm'])
print()

# test exoport
analyze.export(stm, 'data', 'fits', ttl, exfmto='svg', exfmtf='svg')

input("OK")
