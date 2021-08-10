#!/usr/bin/env python
# coding: utf-8

"""
Tool for loading data stored in simulation output files.
"""

from . import *

HL = 8 # number of lines in the header

HQ = { # header quantities and corresponding lines
    'n': 0,
    's': 2,
    'g': 3,
    'z': 4,
    'b': 5,
    'C': 6,
    'a': 7,
}

AVGF = "_averaged.dat" # name extension convention for averaged files

@beartype
def load_file(
    q: List,
    n: str,
    p: str,
) -> Tuple:
    """
    Return the values of the quantities in q from the file n.

    When the requested quantities are provided in the header, only the
    header of the file is loaded.

    Input:
        q: name of the quantities to extract
        n: name of the simulation output file
        p: path to the directory where n can be found

    Output:
        v: values of the quantities in the order requested in q

    The following fields can be extracted:
        'n' (int): number of dislocations
        's' (int): size of the region of interest [nm]
        'g' (Vector): diffraction vector hkl direction
        'z' (Vector): line vector hkl direction
        'b' (Vector): Burgers vector hkl direction
        'C' (Scalar): contrast factor
        'a' (Scalar): latice parameter [nm]
        'J' (int): number of diffraction vector harmonics
        'L' (ScalarList): Fourier variable [nm]
        'cos_AL' (ScalarList): cos of the first harmonic
        'sin_AL' (ScalarList): sin of the first harmonic
        'Cos_<j>AL' (ScalarList): cos of the <j>th harmonic
        'Sin_<j>AL' (ScalarList): sin of the <j>th harmonic
        'err_Cos' (ScalarList): cos error of the first harmonic
        'err_Sin' (ScalarList): sin error of the first harmonic
        'err_Cos_<j>AL' (ScalarList): cos error of the <j>th harmonic
        'err_Sin_<j>AL' (ScalarList): sin error of the <j>th harmonic
        '<eps^2>' (ScalarList): mean square strain
        'bad_points' (ScalarList): number of incorrect random points
    """
    with open(p+n, "r") as f:
        hd = [f.readline().strip('\n') for i in range(HL)] # header data
        tq = f.readline().split() # table quantities
        v, all, i = [], False, 0
        while i<len(q) and not all:
            if q[i] in tq:
                all = True
                td = [l.split() for l in f.readlines()] # table data
            i += 1
    for n in q:
        if n in HQ:
            nv = np.array([eval(v) for v in hd[HQ[n]].split('#')[0].split()])
            if len(nv) == 1:
                nv = nv[0]
            v.append(nv)
        elif n in tq:
            j = tq.index(n)
            v.append(np.array([eval(td[i][j]) for i in range(len(td))]))
        elif n == 'J':
            v.append((len(tq)-3)//4)
        else:
            raise ValueError("unknown quantity: "+n)
    return tuple(v)

@beartype
def load_directory(
    q: List,
    n: str,
    p: str,
) -> Tuple:
    """
    Average and return the values of q over the files in directory n.

    The quantity values of each file in the directory are loaded with
    the load_file function and returned averaged. The same quantities
    as in function load_file can be extracted.

    Input:
        q: name of the quantities to extract and average
        n: name of the simulation output directory
        p: path to the directory where n can be found

    Ouput:
        v: averaged values of the quantities in the order requested
    """
    d = os.listdir(p+n)
    dv = [load_file(q, f, p+n+"/") for f in d]
    v = [sum([dv[i][j] for i in range(len(d))])/len(d) for j in range(len(q))]
    return tuple(v)

@beartype
def average_file(
    s: str,
    i: str,
    o: str,
) -> None:
    """
    Export an averaged simulation output file from a directory.

    Input:
        s: name of the simulation output directory
        i: path where the sample directory can be found
        o: path where the averaged file will be exported
    """
    # recover the header in a random file
    n = os.listdir(i+s)[0]
    with open(i+s+"/"+n, 'r') as f:
        h = [f.readline() for i in range(HL)]
        q = f.readline().split()
    # averaged values
    t = load_directory(q, s, i)
    h[0] = str(load_directory(['n'], s, i)[0])+" #"+h[0].split('#')[1]
    # export file
    fmt = "%5.1f "+" ".join(["%10.7f" for i in range(len(q)-2)])+" %1.1f"
    with open(o+s+AVGF, "w+") as f:
        for l in h:
            f.write(l)
        f.write(" ".join(q)+"\n")
        np.savetxt(f, np.transpose(t), fmt=fmt)

@beartype
def load(
    q: List,
    n: str,
    p: str,
) -> Tuple:
    """
    Return the values of q from the simulation output n.

    When the results of a sample are requested, an average file is
    created if it does not exist. The same quantities as in function
    load_file can be extracted.

    Input:
        q: name of the quantities to extract
        n: name of the simulation output (file or directory)
        p: path where n can be found

    Output:
        v: averaged values of the quantities in the order requested
    """
    if os.path.isfile(p+n):
        return load_file(q, n, p)
    elif os.path.isdir(p+n) or os.path.exists(p+n+AVGF):
        if not os.path.exists(p+n+AVGF):
            average_file(n, p, p)
        return load_file(q, n+AVGF, p)
    else:
        raise ValueError('nothing found at specified path: '+p+n)
