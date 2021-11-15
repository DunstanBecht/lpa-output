#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module collect.
"""

from lpa.output import collect

outfil = 'rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_PBC1_S0_output'
outdir = '10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output'
qtynam = ('v', 'nd', 's', 'd')

# load file
print(collect.load_file(qtynam, outfil, impdir='data'))
print()

# load directory
print(collect.load_directory(qtynam, outdir, impdir='data'))
print()

# versatile load
print(collect.load(qtynam, outfil+'.dat', impdir='data'))
print(collect.load(qtynam, outdir, impdir='data/'))
print()

input("OK")
