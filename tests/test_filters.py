#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module filters.
"""

from lpa.output import filters
import numpy as np
import matplotlib.pyplot as plt

"""
The following lines consist of a visual demonstration of filter 1.
"""
for a in [
    np.array([10, 9.5, 8.5, 6.5, 2.5, 0, -2, -1])/10,
    np.array([10, 9.5, 8.5, 6.5, 2.5, 3, -2, -1])/10,
]:
    l = np.arange(len(a))+1
    i1 = filters.F1(a)
    plt.plot(l, a, "o-", label="all")
    plt.plot(l[:i1], a[:i1], ".-", label="filtered")
    plt.title("F1")
    plt.legend()
    plt.grid()
    plt.show()

"""
The following lines consist of a visual demonstration of filter 2.
"""
for y in [
    np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]),
    np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.2, 6.8, 8.6, 11.2]),
]:
    x = np.arange(len(y))+1
    i2 = filters.F2_xy(y, x)
    plt.plot(x, y, "o-", label="all")
    plt.plot(x[:i2], y[:i2], ".-", label="filtered")
    plt.title("F2")
    plt.legend()
    plt.grid()
    plt.show()

input("OK")
