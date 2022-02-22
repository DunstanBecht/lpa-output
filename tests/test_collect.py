#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module collect.
"""

from lpa.output import collect

outfil = 'rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_PBC1_S0_output'
outdir = '10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output'
qtynam = ('v', 'nd', 's', 'd')

"""
The following line loads the values of the requested quantities from a
file corresponding to a simulation output.
"""
print(collect.load_file(qtynam, outfil, impdir='data'))
print()

"""
The following line generates a file that corresponds to the merging of
the output data of the different distributions of the same sample.
"""
collect.average(outdir, impdir='data')

"""
The following line loads the values of the requested quantities from a
folder corresponding to the simulation output of a distribution sample.
"""
print(collect.load_directory(qtynam, outdir, impdir='data'))
print()

"""
The following line loads the values of the requested quantities from
either a folder or a file.
"""
print(collect.load(qtynam, outfil+'.dat', impdir='data'))
print(collect.load(qtynam, outdir, impdir='data/'))
print()

input("OK")
