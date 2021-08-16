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

@beartype
def load_file(
    q: List,
    imstm: str,
    imdir: str = '',
    imfmt: str = 'dat',
) -> Tuple:
    """
    Return the values of the quantities in q from an output file.

    When the requested quantities are provided in the header, only the
    header of the file is loaded.

    Input:
        q: name of the quantities to extract
        imstm: stem of the file to import
        imdir: import directory
        imfmt: import format

    Output:
        v: values of the quantities in the order requested in q

    The following quantities can be loaded:
        'n' (int): number of dislocations
        's' (int): size of the region of interest [nm]
        'g' (Vector): diffraction vector direction (hkl)
        'z' (Vector): direction of 'l' (line vector) [uvw]
        'b' (Vector): Burgers vector direction [uvw]
        'C' (Scalar): contrast factor [1]
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
        '<eps^2>' (ScalarList): mean square strain [1]
        'bad_points' (ScalarList): number of incorrect random points
    """
    with open(os.path.join(imdir, imstm+'.'+imfmt), "r") as f: # load the file
        hv = [f.readline().strip("\n") for i in range(HL)] # header values
        tq = f.readline().split() # table quantities
        v = [] # values of the quantities in q
        imtab = False # the coefficient table has been imported
        i = 0 # index on q
        while i<len(q) and not imtab: # stops when an element is in the table
            if q[i] in tq: # the quantity must be imported from the table
                td = [l.split() for l in f.readlines()] # load table data
                imtab = True # inform that the table has been imported
            i += 1 # go to next requested quantity
    for n in q:
        if n in HQ: # the quantity is in the header
            nv = np.array([eval(v) for v in hv[HQ[n]].split('#')[0].split()])
            if len(nv) == 1: # if the quantity is a scalar
                nv = nv[0] # return a scalar
            v.append(nv)
        elif n in tq: # the quantity is in the table
            j = tq.index(n)
            v.append(np.array([eval(td[i][j]) for i in range(len(td))]))
        elif n == 'J': # the quantity is the number of harmonics
            v.append((len(tq)-3)//4)
        else:
            raise ValueError("unknown quantity: "+n)
    return tuple(v)

@beartype
def load_directory(
    q: List,
    imstm: str,
    imdir: str = '',
) -> Tuple:
    """
    Average and return the values of q over the files of a directory.

    The quantity values of each file in the directory are loaded with
    the load_file function and returned averaged. The same quantities
    as in function load_file can be extracted.

    Input:
        q: name of the quantities to extract and average
        imstm: stem of the directory to import
        imdir: import directory

    Ouput:
        v: averaged values of the quantities in the order requested
    """
    dv = [] # values for each distribution file
    imdir_stm = os.path.join(imdir, imstm)
    stm_fmt = [os.path.splitext(e) for e in os.listdir(imdir_stm)] # files
    for stm, fmt in stm_fmt:
        dv.append(load_file(q, stm, imdir_stm, fmt[1:])) # load
    v = [] # averaged values
    for j in range(len(q)):
        v.append(sum([dv[i][j] for i in range(len(stm_fmt))])/len(stm_fmt))
    return tuple(v)

@beartype
def average_file(
    imstm: str,
    imdir: str = '',
    exdir: Optional[str] = None,
    exstm: Optional[str] = None,
    exfmt: str = 'dat',
) -> None:
    """
    Export an averaged simulation output file from a directory.

    Input:
        imstm: stem of the directory to average
        imdir: import directory
        exdir: export directory
        exstm: stem of the averaged file
        exfmt: format of the averaged file
    """
    if exdir is None:
        exdir = imdir # default export directory
    if exstm is None:
        exstm = imstm # default export stem
    # recover the header in a random file
    imdir_stm = os.path.join(imdir, imstm)
    n = os.listdir(imdir_stm)[0] # select a file
    with open(os.path.join(imdir_stm, n), 'r') as f: # load the file
        hv = [f.readline() for i in range(HL)] # header values
        tq = f.readline().split() # table quantities
    # average values
    m = load_directory(['n'], imstm, imdir)[0] # averaged number of disl.
    hv[0] = str(m)+" #"+hv[0].split('#')[1] # update the header information
    tv = load_directory(tq, imstm, imdir) # averaged table values
    # export file
    fmt = "%5.1f "+" ".join(["%10.7f" for i in range(len(tq)-2)])+" %1.1f"
    with open(os.path.join(exdir, exstm+'.'+exfmt), "w+") as f:
        for l in hv:
            f.write(l) # write heder lines
        f.write(" ".join(tq)+"\n") # write table
        np.savetxt(f, np.transpose(tv), fmt=fmt)

@beartype
def load(
    q: List,
    imstm: str,
    imdir: str = '',
    imfmt: str = 'dat',
) -> Tuple:
    """
    Return the values of q from a simulation output.

    When the results of a sample are requested, an average file is
    created if it does not exist. The same quantities as in function
    load_file can be extracted.

    Input:
        q: name of the quantities to extract
        imstm: stem of the output to load
        imdir: import directory
        imfmt: import format in case of file output

    Output:
        v: averaged values of the quantities in the order requested
    """
    imdir_stm = os.path.join(imdir, imstm)
    if os.path.isfile(imdir_stm+'.'+imfmt): # load the averaged file
        return load_file(q, imstm, imdir, imfmt)
    elif os.path.isdir(imdir_stm): # average the directory and load file
        average_file(imstm, imdir)
        return load_file(q, imstm, imdir, imfmt)
    else:
        raise ValueError('nothing found at specified path: '+imdir_stm)
