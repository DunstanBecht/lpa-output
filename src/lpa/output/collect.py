#!/usr/bin/env python
# coding: utf-8

"""
Tool for loading data stored in simulation output files.
"""

from . import *

HL = 11 # number of lines in the header

HQ = { # header quantities and corresponding lines
    'v':   0, # lpa-xrd version
    'd':   1, # dislocation density
    'z':   2, # direction of 'l' (line vector)
    'b':   3, # Burgers vector direction
    'g':   4, # diffraction vector direction
    'C':   5, # contrast coefficient
    'a':   6, # cell parameter
    's':   7, # size of the region of interest
    'nu':  8, # Poisson's number
    'nd':  9, # number of dislocations in the input file
    'np': 10, # number of random points
}

@beartype
def load_file(
    qtynam: Union[list, tuple],
    impstm: str,
    **kwargs,
) -> tuple:
    """
    Return the values of the quantities in qtynam from the output file.

    When the requested quantities are provided in the header, only the
    header of the file is loaded.

    Input:
        qtynam (list|tuple): name of the quantities to extract from the file
        impstm (str): stem of the file to import
      **impdir (str): import directory (default: '')
      **impfmt (str): import format (default: 'dat')

    Output:
        qtyval: quantity values in the order of qtynam

    The following variables can be loaded:
        'stm' (str): output stem
        'v' (str): lpa-xrd version
        'd' (Scalar): dislocation density [m^-2]
        'z' (Vector): direction of 'l' (line vector) [uvw]
        'b' (Vector): Burgers vector direction [uvw]
        'g' (Vector): diffraction vector direction (hkl)
        'C' (Scalar): contrast factor [1]
        'a' (Scalar): lattice parameter [nm]
        's' (int): size of the region of interest [nm]
        'nu' (Scalar): Poisson's number [1]
        'nd' (int): number of dislocations in the input file
        'np' (int): number of random points
        'A' (ScalarListList): Fourier transform for each harmonic
        'L' (ScalarList): Fourier variable [nm]
        'cos<j>' (ScalarList): with <j> the harmonic [1]
        'sin<j>' (ScalarList): with <j> the harmonic [1]
        'err_cos<j>' (ScalarList): with <j> the harmonic [1]
        'err_sin<j>' (ScalarList): with <j> the harmonic [1]
        '<eps^2>' (ScalarList): mean square strain [1]
        'bad_points' (ScalarList): number of incorrect random points
    """
    # optional parameters
    impdir = getkwa('impdir', kwargs, str, '')
    impfmt = getkwa('impfmt', kwargs, str, 'dat')
    endkwa(kwargs)
    # load file
    with open(os.path.join(impdir, impstm+'.'+impfmt), "r") as f:
        hv = [f.readline().strip("\n") for i in range(HL)] # header values
        tq = f.readline().split()[1:] # table quantities
        imptab = False # the coefficient table has been imported
        i = 0 # index on qtynam
        while i<len(qtynam) and not imptab: # stops when element is in table
            if qtynam[i]=='A' or qtynam[i] in tq: # the quantity is in table
                td = [l.split() for l in f.readlines()] # load table data
                imptab = True # inform that the table has been imported
            i += 1 # go to next requested quantity
    qtyval = [] # values of the quantities in qtynam
    # sort values
    def aux(nam: str):
        if nam in HQ: # the quantity is in the header
            v = np.array([eval(v) if nam!='v' else v
                          for v in hv[HQ[nam]].split('#')[0].split()])
            if len(v) == 1: # if the quantity is a scalar
                v = v[0] # return a scalar
            return v
        elif nam in tq: # the quantity is in the table
            j = tq.index(nam)
            return np.array([eval(td[i][j]) for i in range(len(td))])
        elif nam == 'A':
            a = (
                np.array([
                    np.array([eval(td[i][4*h+1]) for i in range(len(td))])
                +1j*np.array([eval(td[i][4*h+3]) for i in range(len(td))])
                for h in range((len(tq)-3)//4)])
            )
            return a
        elif nam == 'stm':
            return impstm
        else:
            raise ValueError(f"unknown quantity: {nam}")
    for nam in qtynam:
        qtyval.append(aux(nam))
    return tuple(qtyval)

@beartype
def load_directory(
    qtynam: Union[list, tuple],
    impstm: str,
    **kwargs,
) -> tuple:
    """
    Average and return the values of qtynam over multiple files.

    The quantity values of each file in the directory are loaded with
    the load_file function and returned averaged. The same quantities
    as in function load_file can be extracted.

    Input:
        qtynam (list|tuple): name of the quantities to extract and average
        impstm (str): stem of the directory to import
      **impdir (str): import directory (default: see load_file)

    Ouput:
        qtyval: averaged quantity values in the order of qtynam
    """
    # optional parameters
    impdir = getkwa('impdir', kwargs, str, '')
    endkwa(kwargs)
    # load files
    dv = [] # values for each distribution file
    dir_stm = os.path.join(impdir, impstm) # files directory
    stm_fmt = [os.path.splitext(e) for e in os.listdir(dir_stm)] # files
    for stm, fmt in stm_fmt:
        val = load_file(qtynam, stm, impdir=dir_stm, impfmt=fmt[1:]) # values
        dv.append(val) # store data
    qtyval = [] # averaged values
    for j in range(len(qtynam)):
        if qtynam[j]=='v':
            qtyval.append(dv[0][j])
        elif qtynam[j]=='stm':
            qtyval.append(impstm)
        else:
            lstval = [dv[i][j] for i in range(len(stm_fmt))] # to average
            qtyval.append(sum(lstval)/len(lstval)) # store mean
    return tuple(qtyval)

@beartype
def average(
    impnam: str,
    **kwargs: str,
) -> None:
    """
    Produce a file resulting from the averaging of a sample.

    Input:
        impnam (str): name of the output directory
      **impdir (str): import directory (default: '')
      **expdir (str): export directory (default: impdir)
      **fmtout (str): data export format (default: 'dat')
    """
    impdir = getkwa('impdir', kwargs, str, '')
    expdir = getkwa('expdir', kwargs, str, impdir)
    fmtout = getkwa('fmtout', kwargs, str, 'dat')
    endkwa(kwargs)
    smpdir = os.path.join(impdir, impnam)
    avgfil = os.path.join(expdir, impnam)
    # retrive data
    with open(os.path.join(smpdir, os.listdir(smpdir)[0]), 'r') as f:
        hdr = [f.readline() for i in range(12)]
    clm = hdr[-1].split()[1:]
    qtynam = ['d']+clm
    res = load_directory(qtynam, impnam, impdir=impdir)
    dst = res[0]
    tab = list(res[1:])
    j = (len(clm)-3)//4
    n = len(os.listdir(smpdir))
    # edit data
    hdr[1] = f"{dst:8.2e} #" + hdr[1].split('#')[1]
    for h in range(j):
        for i in [1, 3]:
            tab[1+4*h+i] /= np.sqrt(n)
    # write data
    fmt = ['%6.1f'] + ['%10.7f']*(4*j+1) + ['%10d']
    with open(f"{avgfil}.{fmtout}", 'w') as f:
        f.write("".join(hdr))
        np.savetxt(f, np.array(tab).T, fmt=fmt)

@beartype
def load(
    qtynam: Union[list, tuple],
    impnam: str,
    **kwargs,
) -> tuple:
    """
    Return the values of qtynam from a simulation output file or dir.

    Input:
        qtynam (list|tuple): name of the quantities to extract
        impnam (str): name of the directory or output file to load
      **impdir (str): import directory (default: see load_file)

    Output:
        qtyval: quantity values in the order of qtynam
    """
    # optional parameters
    impdir = getkwa('impdir', kwargs, str, '')
    endkwa(kwargs)
    # load output
    outpth = os.path.join(impdir, impnam) # path to the output
    if os.path.isfile(outpth): # load the output file
        stm, fmt = os.path.splitext(impnam)
        return load_file(qtynam, stm, impdir=impdir, impfmt=fmt[1:])
    elif os.path.isdir(outpth): # load and average the output directory
        return load_directory(qtynam, impnam, impdir=impdir)
    else:
        raise ValueError(f"nothing found at specified path: {outpth}")
