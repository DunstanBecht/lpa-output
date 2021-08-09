#!/usr/bin/env python
# coding: utf-8

"""
Module to facilitate the execution of simulations.
"""

from . import *

@beartype
def make(
    executer=os.system
) -> None:
    """
    Compile the simulation program.

    Input:
        executer:
    """
    print(executer("module load cuda/10.1; cd xray-diffraction; make"))

@beartype
def run(
    s: str,
    i: str = "input/",
    o: str = "output/",
    h: int = 1,
    b: int = 200,
    r: int = 1000,
    f: int = 35,
    executer=os.system,
) -> None:
    """
    Run the simulation on a sample of distributions.

    Input:
        s: directory name of the sample of distributions
        i: input directory
        o: output directory
        h: hardware to use (1 for gpu / 0 for cpu)
        b: block size
        r: block repetitions (r*b gives the number of random points)
        f: number of Fourier coefficients
        executer:
    """
    n = r*b # number of random points
    if not os.path.exists(o+s):
        os.mkdir(o+s)
    ipath = i+s+"/"
    opath = o+s+"/"
    for e in os.listdir(i+s):
        args = " ".join(str(a) for a in [h, b, ipath+e, n, f, opath+e])
        print("- "+e+" ("+args+")")
        executer("cd xray-diffraction; ./a.out "+args+" >& runs/"+e)
