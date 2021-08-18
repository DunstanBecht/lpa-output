#!/usr/bin/env python
# coding: utf-8

"""
Tools for the analysis of the simulation output.
"""

import matplotlib
import matplotlib.pyplot as plt
import scipy.optimize
from . import *
from . import collect
from . import models

matplotlib.use("Agg")

@beartype
def A(
    imstm: str,
    imdir: str,
    j: ScalarList,
) -> ScalarListList:
    """
    Return the averaged Fourier amplitude for the harmonics j.

    Input:
        imstm: stem of the simulation output
        j: harmonics of diffraction vector
        imdir: import directory

    Output:
        a: averaged Fourier amplitude for each harmonic in j

    Output example:
        a = [
            np.array([A_1(L_1), A_1(L_2), A_1(L_3), ...]),
            np.array([A_2(L_1), A_2(L_2), A_2(L_3), ...]),
            ...
        ]
    """
    q = [] # quantity names
    for h in j:
        if h == 1:
            q += ['cos_AL', 'sin_AL']
        else:
            q += ['Cos_'+str(h)+'AL', 'Sin_'+str(h)+'AL']
    m = collect.load(q, imstm, imdir) # load coefficidnt table
    return [np.sqrt(m[2*h-2]**2+m[2*h-1]**2) for h in j]

@beartype
def lowpass1(
    a: ScalarList,
) -> int:
    """
    Return the index from which the noise begins.

    Input:
        a: Fourier amplitudes for a given harmonic

    Output:
        i1: index at which the noise starts
    """
    i = 1
    while i<len(a)-1 and a[i-1]>a[i]:
        i += 1
    return min(len(a)-1, i+4)

@beartype
def lowpass2(
    a: ScalarList,
    l: ScalarList,
) -> int:
    """
    Return the index from which ln(A(L))/L^2 is no longer linear.

    Input:
        l: Fourier variable
        a: Fourier amplitudes for a given harmonic

    Output:
        i2: index that marks the end of the linear part
    """
    errors = []
    y = np.log(a)/l**2
    x = np.stack((np.log(l), np.ones(len(l)))).T
    #plt.plot(x[:,0], y)
    #plt.show()
    for i in range(4, len(l)+1):
        errors.append(np.linalg.lstsq(x[:i], y[:i], rcond=None)[1][0]/i)
    #print(errors)
    #plt.plot(l[3:], errors)
    #plt.show()
    m = max(errors)
    for i in range(len(errors)):
        if errors[i]/m>0.001:
            return i + 3
    return 0

@beartype
def common(
    imstm: str,
    imdir: str,
    j: Optional[ScalarList] = None,
) -> dict:
    """
    Return a set of common quantities related to a simulation.

    Input:
        imstm: stem of the simulation output
        imdir: import directory
        j: harmonics studied

    Output:
        c: dictionary containing the quantities common to a simulation

    The quantities contained in c are the following:
        'name': name of the simulation output
        'g' (Vector): diffraction vector [nm^-1]
        'z' (Vector): line vector direction [uvw]
        'b' (Vector): Burgers vector [nm]
        'C' (Scalar): contrast factor [1]
        'a' (Scalar): latice parameter [nm]
        'J' (int): number of diffraction vector harmonics
        'L' (ScalarList): Fourier variable [nm]
        'j' (ScalarList): harmonics studied
        'jg' (VectorList): j*g [nm^-1]
        'jgb' (ScalarList): (j*g).b [1]
        'b2' (Scalar): |b|^2 [nm^2]
        'jg2' (Scalar): |j*g|^2 [nm^-2]
        'A' (ScalarListList): Fourier amplitudes for harmonics in j [1]
        'i1' (List): index at which the noise starts for each harmonic
        'i2' (List): index of linear zone end for each harmonic
        'index' (dict): reverse indexing for the harmonics
    """
    c = {}
    q = ['g', 'z', 'b', 'C', 'a', 'J', 'L']
    for k, v in zip(q, collect.load(q, imstm, imdir)):
        c[k] = v
    c['g'] = c['g']/c['a'] # diffraction vector with correct norm
    if j is None:
        c['j'] = np.arange(c['J']) + 1
    else:
        c['j'] = j
    c['jg'] = c['j'].reshape((len(c['j']), 1))*c['g']
    c['jgb'] = np.sum(c['b']*c['jg'], axis=1)
    c['b2'] = np.sum(c['b']**2)
    c['jg2'] = np.sum(c['jg']**2, axis=1)
    c['A'] = A(imstm, imdir, c['j'])
    c['i1'] = [lowpass1(a) for a in c['A']]
    c['i2'] = [lowpass2(a, c['L']) for a in c['A']]
    c['index'] = {}
    for i in range(len(c['j'])):
        c['index'][c['j'][i]] = i
    return c

@beartype
def fit(
    m: Callable,
    c: dict,
    f: str,
    d: Scalar = 5e14*1e-18,
    r: Scalar = 200,
) -> dict:
    """
    Return information on the fits made of model m on simulation c.

    A fit is calculated for each harmonic and each interval of L.

    Input:
        m: model function
        c: dictionary containing the quantities common to a simulation
        f: lowpass filter name
        d: initial density [nm^-2]
        r: initial outer cut-off radius [nm]

    Output:
        f: dictionary containing the information on the fits

    f contains the following fields:
        'name' (str): name of the model
        'j' (ScalarList): harmonics
        'L' (ScalarList): maximum values ​​of L [nm]
        'd' (ScalarList): optimal density for the fit [nm^-2]
        'r' (ScalarList): optimal outer cut-off radius for the fit [nm]
        'e' (ScalarList): optimal fit errors
        'm' (Callable): model function
    """
    J, L, D, R, E, = [], [], [], [], [] # fits information
    for i_j in range(len(c['j'])):
        for i_L in range(3, c[f][i_j]):
            j = c['j'][i_j] # harmonic
            l = c['L'][:i_L] # maximum value ​​of L
            a = c['A'][i_j][:i_L] # Fourier amplitudes
            def error(p)-> Scalar:
                return np.sum((a-m(*p, c, int(j), l))**2/(i_L-2))
            print('')
            p = scipy.optimize.fmin(error, (d, r), ftol=1e-10)
            print('')
            J.append(j) # add harmonic
            L.append(l[-1]) # add maximum value ​​of L
            D.append(p[0]) # add density
            R.append(p[1]) # add outer cut-off radius
            E.append(error(p)) # add error
    f = {
        'name': m.__name__,
        'j': np.array(J),
        'L': np.array(L),
        'd': np.array(D),
        'r': np.array(R),
        'e': np.array(E),
        'm': m,
    }
    return f

@beartype
def plot(
    c: dict,
    exstm: str,
    exdir: str,
    exfmt: str = 'pdf',
    title: str = '',
    j: Optional[ScalarList] = None,
    f: dict = None,
    L: Optional[ScalarList] = None,
) -> None:
    """
    Export a figure with two representations of A(L) and the fits.

    Input:
        c: common parameters dictionary
        exstm: export stem
        exdir: export directory
        exfmt: export format
        title: title of the figure
        j: restriction of harmonics to be displayed
        f: fits information
        L: restriction of maximum L values
    """
    if j is None: # no restriction of the harmonics to be displayed is applied
        j = c['j'] # display all the available harmonics
    data = [] # list of curves to display
    for h in j:
        i_j = c['index'][h]
        i_L = c['i1'][i_j]
        data.append({
            'L': c['L'][:i_L], # x variable
            'A': c['A'][i_j][:i_L], # y variable
            'm': '.', # marker
            'n': "A_{"+str(h)+"}(L)", # label
            'c': 'C'+str(h-1), # color
            'l': '', # extra information
            })
    if not f is None: # display fits
        mask = np.in1d(f['j'], j) # restrict to the requested harmonics
        if not L is None: # restrict the Fourier variable range
            mask = mask & np.in1d(f['L'], L) # # restrict to the requested L
        maskj = f['j'][mask] # apply mask on j
        maskL = f['L'][mask] # apply mask on L
        maskd = f['d'][mask] # apply mask on d
        maskr = f['r'][mask] # apply mask on r
        maske = f['e'][mask] # apply mask on e
        min = c['L'][0] # minimum value of L (common to all fits)
        m = f['m'] # model function
        for i in range(len(maskj)): # through remaining fits
            max = maskL[i] # maximum value of L
            d, r = maskd[i], maskr[i] # density and outer cut-off radius
            v = r" \times 10^{".join(format(maske[i], '1.1e').split('e'))+"}"
            le = r"$ \sigma^2 ="+v+" $" # fit error
            ll = r"$ L \leq "+str(max)+r" $" # fit range
            h = maskj[i] # harmonic
            l = np.linspace(min, max, 40) # range to plot the fit
            data.append({
                'L': l, # x variable
                'A': m(d, r, c, int(h), l), # y variable
                'm': '-', # marker
                'n': "A^"+f['name'][0]+"_{"+str(h)+r"}(L)", #  label
                'c': 'C'+str(h-1), # color
                'l': ", "+le+", "+ll, # extra information
            })
    # fig
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.1)
    fig.suptitle(title)
    # ax1: A(L) as a function of L (log scale)
    for d in data:
        ax1.plot(
            d['L'],
            d['A'],
            d['m'],
            label=r"$"+d['n']+r"$"+d['l'],
            color=d['c'],
        )
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$L \ (nm)$")
    ax1.grid()
    ax1.legend()
    # ax2: ln(A(L)) (log scale) as a function of L^2
    for d in data:
        ax2.plot(
            d['L'],
            np.log(d['A'])/d['L']**2,
            d['m'],
            label=r"$\frac{\ln("+d['n']+r")}{L^2}$"+d['l'],
            color=d['c'],
        )
    ax2.set_xscale("log")
    ax2.set_xlabel(r"$L \ (nm)$")
    ax2.grid()
    ax2.legend()
    # export
    plt.savefig(exdir+exstm+'.'+exfmt, format=exfmt)
    fig.clf()

@beartype
def export_model(
    m: dict,
    c: dict,
    exdir: str,
    exstm: str,
    exfmtd: str = 'csv',
    exfmtp: str = 'png',
    d: Scalar = 5e14*1e-18,
    r: Scalar = 200,
) -> None:
    """
    Export the information on fits and the corresponding figures.

    Input:
        m: dictionary containing information about the model
        c: common quantities
        exdir: export directory
        exstm: export stem
        exfmtd: data export format
        exfmtp: plot export format
        d: initial density [nm^-2]
        r: initial outer cut-off radius [nm]
    """
    f = fit(m['function'], c, m['filter']) # perform the fits
    # export a figure for each fit
    exdir_mod = exdir+"/"+f['name']+"/" # model export directory
    if not os.path.exists(exdir_mod):
        os.mkdir(exdir_mod)
    for i in range(len(f['j'])):
        fit_stm = ("j"
            + format(f['j'][i], '01.0f')
            + "_"
            + format(f['L'][i], '03.0f')
            + "nm"
        )
        plot(
            c,
            fit_stm,
            exdir_mod,
            exfmt=exfmtp,
            f=f,
            j=np.array([f['j'][i]]),
            L=np.array([f['L'][i]]),
        )
    # export fits data
    fields = (
        'harmonic of g',
        'fit L max [nm]',
        'rho [m-2]',
        'Re [nm]',
        'error',
    )
    values = (f['j'], f['L'], f['d']*1e18, f['r'], f['e'],)
    sep = ";"
    fmt = ['%1.0f', '%3.1f', '%1e', '%1e', '%1e']
    with open(exdir+exstm+"_"+f['name']+"."+exfmtd, "w") as f:
        f.write(sep.join(fields)+'\n')
        np.savetxt(f, np.transpose(values), fmt=sep.join(fmt))

@beartype
def export(
    imstm: str,
    imdir: str = "",
    exdir: str = "",
    title: Optional[str] = None,
    d: Scalar = 5e14*1e-18,
    r: Scalar = 200,
) -> None:
    """
    Perform an analysis with the available models.

    Input:
        imstm: stem the simulation output (file or directory)
        imdir: input directory
        exdir: export directory
        title: figure title
        d: initial density [nm^-2]
        r: initial outer cut-off radius [nm]
    """
    if imdir!="" and imdir[-1]!="/":
        imdir += "/"
    if exdir!="" and exdir[-1]!="/":
        exdir += "/"
    if title is None:
        title = " ".join(imstm.split("_"))
    # directory
    exdir_stm = exdir+imstm+"/"
    if not os.path.exists(exdir_stm):
        os.mkdir(exdir_stm)
    # general
    c = common(imstm, imdir, np.array([1, 2]))
    plot(c, imstm, exdir_stm, title=title)
    # models
    analyzed = [
        {'function': models.Groma, 'filter': 'i2',},
        {'function': models.Kamminga, 'filter': 'i2',},
        {'function': models.Wilkens, 'filter': 'i1',},
    ]
    for m in analyzed:
        export_model(m, c, exdir_stm, imstm, d=d, r=r)
