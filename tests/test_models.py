#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module analyze.
"""

from lpa.output import analyze, models
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('TkAgg')

outdat = {
    'j': np.array((1, 2)),
    'g': np.array((2, 0, 0)),
    'b': np.array((0.202470, 0.202470, 0.000000)),
    'z': np.array((1, 1, 0)),
    'L': np.linspace(1, 300, 100),
    'C': 0.250000,
}
vertij = outdat['j'].reshape((len(outdat['j']), 1))
outdat['jg'] = vertij*outdat['g']
outdat['jgb'] = np.sum(outdat['b']*outdat['jg'], axis=1)
outdat['b2'] = np.sum(outdat['b']**2)
outdat['jg2'] = np.sum(outdat['jg']**2, axis=1)
outdat['index'] = {outdat['j'][i]: i for i in range(len(outdat['j']))}

h = 1

# GUW1
p = np.array([5e-05, 4500, 1.5, 1000])
plt.plot(outdat['L'], models.GUW1(p, outdat, h, outdat['L']), label='GUW1')

# GUW2
p = np.array([5e-05, 4500])
plt.plot(outdat['L'], models.GUW2(p, outdat, h, outdat['L']), label='GUW2')

# W1
p = np.array([5e-05, 500])
plt.plot(outdat['L'], models.W1(p, outdat, h, outdat['L']), label='W1')

# W2
p = np.array([5e-05, 500])
plt.plot(outdat['L'], models.W2(p, outdat, h, outdat['L']), label='W2')

plt.legend()
plt.show()
