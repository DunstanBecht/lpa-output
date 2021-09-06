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
from . import filters

matplotlib.use("Agg") # to export plots with no bitmap allocation errors

@beartype
def output_data(
    imstm: str,
    imdir: str,
) -> dict:
    """
    Return a set of quantities related to a simulation output.

    Input:
        imstm: stem of the simulation output
        imdir: import directory

    Output:
        o: output data dictionary

    The quantities contained in o are the following:
        'stm': stem of the simulation output
        'g' (Vector): diffraction vector [nm^-1]
        'z' (Vector): line vector direction [uvw]
        'b' (Vector): Burgers vector [nm]
        'C' (Scalar): contrast factor [1]
        'a' (Scalar): latice parameter [nm]
        'L' (ScalarList): Fourier variable [nm]
        'A' (ScalarListList): Fourier amplitudes for harmonics in j [1]
        'j' (ScalarList): harmonics of diffraction vector
        'jg' (VectorList): j*g [nm^-1]
        'jgb' (ScalarList): (j*g).b [1]
        'b2' (Scalar): |b|^2 [nm^2]
        'jg2' (Scalar): |j*g|^2 [nm^-2]
        'f1' (List): index at which the noise starts for each harmonic
        'f2' (List): index of linear zone end for each harmonic
        'index' (dict): reverse indexing for the harmonics
    """
    o = {'stm': imstm}
    q = ['g', 'z', 'b', 'C', 'a', 'L', 'A'] # quantities to load
    for k, v in zip(q, collect.load(q, imstm, imdir)):
        o[k] = v # store the loaded quantities
    o['g'] = o['g']/o['a'] # correct diffraction vector norm
    o['j'] = np.arange(len(o['A'])) + 1 # harmonics
    # pre-calculated quantities
    o['jg'] = o['j'].reshape((len(o['j']), 1))*o['g'] # j*g [nm^-1]
    o['jgb'] = np.sum(o['b']*o['jg'], axis=1) # (j*g).b [1]
    o['b2'] = np.sum(o['b']**2) # |b|^2 [nm^2]
    o['jg2'] = np.sum(o['jg']**2, axis=1) # |j*g|^2 [nm^-2]
    # filters
    o['f1'] = [filters.f1(a) for a in o['A']] # suppress noise
    o['f2'] = [filters.f2(a, o['L']) for a in o['A']] # keep linear area
    # reverse index directory
    o['index'] = {o['j'][i]: i for i in range(len(o['j']))}
    return o

@beartype
def error(
    p: ScalarList,
    o: dict,
    j: int,
    l: ScalarList,
    a: ScalarList,
    m: ModelFunction,
) -> Scalar:
    """
    Return the error of the fit of m on a.

    Input:
        p: tuple of the density [nm^-2] and the outer cut-off radius [nm]
        o: output data dictionary
        j: harmonic fitted
        l: Fourier variable [nm]
        a: values of the Fourier amplitude
        m: function of the model fitted

    Output:
        e: error of the fit (sigma^2)
    """
    if m in [models.Groma, models.Kamminga]:
        output = np.log(a)/l**2
        model = np.log(m(*p, o, j, l))/l**2
    elif m in [models.Wilkens]:
        output = a
        model = m(*p, o, j, l)
    return np.sum((output-model)**2)/(len(l)-2)

@beartype
def fits_data(
    m: Callable,
    o: dict,
    f: str,
    d: Scalar = 5e14*1e-18,
    r: Scalar = 200,
    b: Tuple = ((5e10*1e-18, 5e18*1e-18), (1, np.inf)),
) -> dict:
    """
    Return information on the fits made of model m on output o.

    A fit is calculated for each harmonic and each interval of L.

    Input:
        m: model function
        o: output data dictionary
        f: lowpass filter name
        d: initial density [nm^-2]
        r: initial outer cut-off radius [nm]
        b: bounds for (d, r)

    Output:
        f: fits data dictionary

    The quantities contained in f are the following:
        'name' (str): name of the model
        'j' (ScalarList): harmonics
        'L' (ScalarList): maximum values ​​of L [nm]
        'd' (ScalarList): optimal density for the fit [nm^-2]
        'r' (ScalarList): optimal outer cut-off radius for the fit [nm]
        'e' (ScalarList): optimal fit errors
        'm' (Callable): model function
    """
    J, L, D, R, E, = [], [], [], [], [] # fits information
    for i_j in range(len(o['j'])): # harmonic index
        i_L = 3 # index of the maximum value of L for the fit
        failed = False # a fit has failed
        while not failed and i_L<=o[f][i_j]:
            j = o['j'][i_j] # harmonic
            l = o['L'][:i_L] # maximum value ​​of L
            a = o['A'][i_j][:i_L] # Fourier amplitudes for harmonic j
            res = scipy.optimize.minimize(
                error,
                np.array((d, r)),
                args=(o, int(j), l, a, m),
                method='Nelder-Mead',
                bounds=b,
                #options={'maxiter': 1e6},
                #tol = 1e-18
            )
            if res.success:
                J.append(j) # store harmonic
                L.append(l[-1]) # store maximum value ​​of L [nm]
                D.append(res.x[0]) # store optimal density [nm^-2]
                R.append(res.x[1]) # store optimal outer cut-off radius [nm]
                E.append(res.fun) # store error
            else:
                print(("/!\ "+res.message+"\n    "
                    + "L="+str(l[-1])+"nm "
                    + "j="+str(j)+" "
                    + "model="+m.__name__
                ))
                failed = True
            i_L += 1
    f = {
        'name': m.__name__,
        'j': np.array(J), # harmonics
        'L': np.array(L), # maximum values of L [nm]
        'd': np.array(D), # optimal densities [nm^-2]
        'r': np.array(R), # optimal outer cut-off radii [nm]
        'e': np.array(E), # errors
        'm': m, # model function
    }
    return f

@beartype
def plot(
    o: dict,
    exstm: str,
    exdir: str,
    exfmt: str = 'pdf',
    title: Optional[str] = None,
    j: Optional[ScalarList] = None,
    f: dict = None,
    L: Optional[ScalarList] = None,
) -> None:
    """
    Export a figure with two representations of A(L) and the fits.

    Input:
        o: output data dictionary
        exstm: export stem
        exdir: export directory
        exfmt: export format
        title: title of the figure
        j: restriction of harmonics to be displayed
        f: fits data dictionary
        L: restriction of maximum L value of the fits to be displayed
    """
    if not title:
        title = exstm.replace("_", " ")
    if j is None: # no restriction of the harmonics to be displayed is applied
        j = o['j'] # display all the available harmonics
    data = [] # list of the information on the curves to display
    for h in j:
        i_j = o['index'][h] # index of the harmonic in c
        i_L = min(o['f1'][i_j]+5, len(o['L'])) # output range to be displayed
        data.append({
            'L': o['L'][:i_L], # x variable
            'A': o['A'][i_j][:i_L], # y variable
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
        Lmin = o['L'][0] # minimum value of L (common to all fits)
        m = f['m'] # model function
        for i in range(len(maskj)): # through remaining fits
            Lmax = maskL[i] # maximum value of L [nm]
            d, r = maskd[i], maskr[i] # density and outer cut-off radius
            v = r" \times 10^{".join(format(maske[i], '1.1e').split('e'))+"}"
            le = r"$ \sigma^2 ="+v+" $" # fit error
            ll = r"$ L \leq "+format(Lmax, '1.1f')+r" $" # fit range
            h = maskj[i] # harmonic
            l = np.linspace(Lmin, Lmax, 40) # range to plot the fit
            data.append({
                'L': l, # x variable
                'A': m(d, r, o, int(h), l), # y variable
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
            d['L'], # x variable
            d['A'], # y variable
            d['m'], # marker
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
            d['L'], # x variable
            np.log(d['A'])/d['L']**2, # y variable
            d['m'], # marker
            label=r"$\frac{\ln("+d['n']+r")}{L^2}$"+d['l'],
            color=d['c'],
        )
    ax2.set_xscale("log")
    ax2.set_xlabel(r"$L \ (nm)$")
    ax2.grid()
    ax2.legend()
    # export
    plt.savefig(os.path.join(exdir, exstm+'.'+exfmt), format=exfmt)
    plt.close('all')

@beartype
def export_model(
    m: dict,
    o: dict,
    exdir: str,
    exstm: str,
    exfmtd: str = 'csv',
    exfmtf: str = 'png',
    title: Optional[str] = None,
    d: Scalar = 5e14*1e-18,
    r: Scalar = 200,
) -> None:
    """
    Export the information on fits and the corresponding figures.

    Input:
        m: dictionary containing information about the model
        o: output data dictionary
        exdir: export directory
        exstm: export stem
        exfmtd: data export format
        exfmtf: fits export format
        title: title of the plots
        d: initial density [nm^-2]
        r: initial outer cut-off radius [nm]
    """
    f = fits_data(m['function'], o, m['filter'], d=d) # get fits information
    # export a figure for each fit
    exdir_mod = exdir+"/"+f['name']+"/" # model export directory
    if not os.path.exists(exdir_mod):
        os.mkdir(exdir_mod)
    for i in range(len(f['e'])): # through the fits
        fit_stm = ("j"
            + format(f['j'][i], '01.0f')
            + "_"
            + format(f['L'][i], '03.0f')
            + "nm"
        )
        plot(
            o, # output data
            fit_stm, # stem of the fit
            exdir_mod, # model export directory
            exfmt=exfmtf, # plot export format
            title=title, # plot title
            f=f, # fits
            j=np.array([f['j'][i]]), # restriction of the harmonic
            L=np.array([f['L'][i]]), # restriction of the maximum value of L
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
    with open(os.path.join(exdir, exstm+"_"+f['name']+"."+exfmtd), "w") as f:
        f.write(sep.join(fields)+'\n')
        np.savetxt(f, np.transpose(values), fmt=sep.join(fmt))

@beartype
def export(
    imstm: str,
    imdir: str = "",
    exdir: str = "",
    title: Optional[str] = None,
    exfmtd: str = 'csv',
    exfmto: str = 'pdf',
    exfmtf: str = 'png',
    d: Scalar = 5e14*1e-18,
    r: Scalar = 200,
    j: ScalarList = np.array([1, 2]),
) -> None:
    """
    Perform an analysis with the available models.

    Input:
        imstm: stem the simulation output (file or directory)
        imdir: input directory
        exdir: export directory
        title: figure title
        exfmtd: data export format
        exfmto: output export format
        exfmtf: fits export format
        d: initial density [nm^-2]
        r: initial outer cut-off radius [nm]
    """
    if title is None:
        title = " ".join(imstm.split("_"))
    # export directory
    exdir_stm = os.path.join(exdir, imstm)
    if not os.path.exists(exdir_stm):
        os.mkdir(exdir_stm)
    # load output data
    o = output_data(imstm, imdir)
    # plot output data
    plot(o, imstm, exdir_stm, title=title, exfmt=exfmto)
    # models
    analyzed = [
        {'function': models.Groma, 'filter': 'f2'},
        {'function': models.Kamminga, 'filter': 'f2'},
        {'function': models.Wilkens, 'filter': 'f1'},
    ]
    # fits
    for m in analyzed:
        export_model(
            m,
            o,
            exdir_stm,
            imstm,
            d=d,
            r=r,
            title=title,
            exfmtd=exfmtd,
            exfmtf=exfmtf,
        )
