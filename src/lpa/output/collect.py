#!/usr/bin/env python
# coding: utf-8

"""
Tool for loading data stored in simulation output files.
"""

from . import *

HL = 8 # number of lines in the header

HQ = { # header quantities and corresponding lines
    'n': 0, # number of dislocations
    's': 2, # size of the region of interest [nm]
    'g': 3, # diffraction vector direction (hkl)
    'z': 4, # direction of 'l' (line vector) [uvw]
    'b': 5, # Burgers vector direction [uvw]
    'C': 6, # contrast factor [1]
    'a': 7, # latice parameter [nm]
}

@beartype
def load_file(
    qtynam: Union[List, Tuple],
    impstm: str,
    **kwargs,
) -> Tuple:
    """
    Return the values of the quantities in qtynam from the output file.

    When the requested quantities are provided in the header, only the
    header of the file is loaded.

    Input:
        qtynam: name of the quantities to extract from the file
        impstm: stem of the file to import
      **impdir: import directory (default: '')
      **impfmt: import format (default: 'dat')

    Output:
        qtyval: quantity values in the order of qtynam

    The following quantities can be loaded:
        'n' (int): number of dislocations
        's' (int): size of the region of interest [nm]
        'g' (Vector): diffraction vector direction (hkl)
        'z' (Vector): direction of 'l' (line vector) [uvw]
        'b' (Vector): Burgers vector direction [uvw]
        'C' (Scalar): contrast factor [1]
        'a' (Scalar): latice parameter [nm]
        'L' (ScalarList): Fourier variable [nm]
        'A' (ScalarListList): Fourier transform for each harmonic
        '<eps^2>' (ScalarList): mean square strain [1]
        'bad_points' (ScalarList): number of incorrect random points
    """
    # optional parameters
    impdir = kwargs.pop('impdir', '') # import directory
    impfmt = kwargs.pop('impfmt', 'dat') # import format
    if len(kwargs)>0:
        raise ValueError("wrong keywords: "+str(kwargs))
    # load file
    with open(os.path.join(impdir, impstm+'.'+impfmt), "r") as f:
        hv = [f.readline().strip("\n") for i in range(HL)] # header values
        tq = f.readline().split() # table quantities
        imptab = False # the coefficient table has been imported
        i = 0 # index on qtynam
        while i<len(qtynam) and not imptab: # stops when element is in table
            if qtynam[i]=='A' or qtynam[i] in tq: # the quantity is in table
                td = [l.split() for l in f.readlines()] # load table data
                imptab = True # inform that the table has been imported
            i += 1 # go to next requested quantity
    qtyval = [] # values of the quantities in qtynam
    # sort values
    for nam in qtynam:
        if nam in HQ: # the quantity is in the header
            v = np.array([eval(v) for v in hv[HQ[nam]].split('#')[0].split()])
            if len(v) == 1: # if the quantity is a scalar
                v = v[0] # return a scalar
            qtyval.append(v)
        elif nam in tq: # the quantity is in the table
            j = tq.index(nam)
            qtyval.append(np.array([eval(td[i][j]) for i in range(len(td))]))
        elif nam == 'A':
            qtyval.append(
                np.array([
                    np.array([eval(td[i][4*h+1]) for i in range(len(td))])**2
                +1j*np.array([eval(td[i][4*h+3]) for i in range(len(td))])**2
                for h in range((len(tq)-3)//4)])
            )
        else:
            raise ValueError("unknown quantity: "+nam)
    return tuple(qtyval)

@beartype
def load_directory(
    qtynam: Union[List, Tuple],
    impstm: str,
    **kwargs,
) -> Tuple:
    """
    Average and return the values of qtynam over multiple files.

    The quantity values of each file in the directory are loaded with
    the load_file function and returned averaged. The same quantities
    as in function load_file can be extracted.

    Input:
        qtynam: name of the quantities to extract and average
        impstm: stem of the directory to import
      **impdir: import directory (default: see load_file)

    Ouput:
        qtyval: averaged quantity values in the order of qtynam
    """
    # optional parameters
    impdir = kwargs.pop('impdir', '') # import directory
    if len(kwargs)>0:
        raise ValueError("wrong keywords: "+str(kwargs))
    # load files
    dv = [] # values for each distribution file
    dir_stm = os.path.join(impdir, impstm) # files directory
    stm_fmt = [os.path.splitext(e) for e in os.listdir(dir_stm)] # files
    for stm, fmt in stm_fmt:
        val = load_file(qtynam, stm, impdir=dir_stm, impfmt=fmt[1:]) # values
        dv.append(val) # store data
    qtyval = [] # averaged values
    for j in range(len(qtynam)):
        lstval = [dv[i][j] for i in range(len(stm_fmt))] # values to average
        qtyval.append(sum(lstval)/len(lstval)) # store mean
    return tuple(qtyval)

@beartype
def load(
    qtynam: Union[List, Tuple],
    impnam: str,
    **kwargs,
) -> Tuple:
    """
    Return the values of qtynam from a simulation output file or dir.

    Input:
        qtynam: name of the quantities to extract
        impnam: name of the directory or output file to load
      **impdir: import directory (default: see load_file)

    Output:
        qtyval: quantity values in the order of qtynam
    """
    # optional parameters
    impdir = kwargs.pop('impdir', '') # import directory
    if len(kwargs)>0:
        raise ValueError("wrong keywords: "+str(kwargs))
    # load output
    outpth = os.path.join(impdir, impnam) # path to the output
    if os.path.isfile(outpth): # load the output file
        stm, fmt = os.path.splitext(impnam)
        return load_file(qtynam, stm, impdir=impdir, impfmt=fmt[1:])
    elif os.path.isdir(outpth): # load and average the output directory
        return load_directory(qtynam, impnam, impdir=impdir)
    else:
        raise ValueError("nothing found at specified path: "+outpth)
