#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module collect.
"""

from lpa.output import collect

output = '10_rho5e14m-2_square_2000nm_RDD_d5e-4nm-2_edge_PBCR2_S0'
q = ['n', 's']

# load file
print(collect.load_file(q, '01', 'data/'+output))
print()

# load directory
print(collect.load_directory(q, output, 'data'))
print()

# versatile load
print(collect.load(q, output, 'data'))
print(collect.load(q, '01.dat', 'data/'+output))
print()

input("OK")
