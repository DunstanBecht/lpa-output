#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module collect.
"""

from lpa.output import collect

output = '10_rho5e14m-2_square_2000nm_RDD_d5e-4nm-2_edge_PBCR2_S0'
qtynam = ('n', 's')

# load file
print(collect.load_file(qtynam, '01', impdir='data/'+output))
print()

# load directory
print(collect.load_directory(qtynam, output, impdir='data'))
print()

# versatile load
print(collect.load(qtynam, output, impdir='data'))
print(collect.load(qtynam, '01.dat', impdir='data/'+output))
print()

input("OK")
