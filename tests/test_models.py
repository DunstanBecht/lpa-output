#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module analyze.
"""

from lpa.output import analyze, models
import numpy as np

stm = '10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output'
outdat = analyze.output_data(stm, impdir='data')
j = 4

# GUW1
p = np.array([4.50138457e-05, 4.59139387e+02, -1.00000000e+00, 1.05716024e+02])
print('GUW1')
print(models.GUW1(p, outdat, j, outdat['L']))
print()

# GUW2
p = np.array([4.50138457e-05, 4.59139387e+02])
print('GUW2')
print(models.GUW2(p, outdat, j, outdat['L']))
print()

# W1
p = np.array([4.50138457e-05, 4.59139387e+02])
print('W1')
print(models.GUW2(p, outdat, j, outdat['L']))
print()

# W2
p = np.array([4.50138457e-05, 4.59139387e+02])
print('W2')
print(models.GUW2(p, outdat, j, outdat['L']))
print()

input("OK")
