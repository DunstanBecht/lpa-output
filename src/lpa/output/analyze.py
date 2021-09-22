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
    impstm: str,
    **kwargs,
) -> dict:
    """
    Return a set of quantities related to a simulation output.

    Input:
        impstm (stm): stem of the simulation output
      **impdir: import directory (default: see collect.load_file)

    Output:
        outdat: output data dictionary

    The quantities contained in outdat are the following:
        'stm' (str): stem of the simulation output
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
    outdat = {'stm': impstm} # output data dictionary
    qtynam = ['g', 'z', 'b', 'C', 'a', 'L', 'A'] # quantities to load
    for key, val in zip(qtynam, collect.load(qtynam, impstm, **kwargs)):
        outdat[key] = val # store the loaded quantities
    outdat['g'] = outdat['g']/outdat['a'] # correct diffraction vector norm
    outdat['j'] = np.arange(len(outdat['A'])) + 1 # harmonics
    # pre-calculated quantities
    vertij = outdat['j'].reshape((len(outdat['j']), 1))
    outdat['jg'] = vertij*outdat['g'] # j*g [nm^-1]
    outdat['jgb'] = np.sum(outdat['b']*outdat['jg'], axis=1) # (j*g).b [1]
    outdat['b2'] = np.sum(outdat['b']**2) # |b|^2 [nm^2]
    outdat['jg2'] = np.sum(outdat['jg']**2, axis=1) # |j*g|^2 [nm^-2]
    # filters
    outdat['f1'] = [filters.f1(a) for a in outdat['A']] # suppress noise
    outdat['f2'] = [filters.f2(a, outdat['L']) for a in outdat['A']] # linear
    # reverse index directory
    outdat['index'] = {outdat['j'][i]: i for i in range(len(outdat['j']))}
    return outdat

@beartype
def standard_error(
    p: ScalarList,
    o: dict,
    j: int,
    l: ScalarList,
    a: ScalarList,
    m: ModelFunction,
    W: bool,
) -> Scalar:
    """
    Return the standard error of the fit of m on a for harmonic j.

    Input:
        p: tuple of the density [nm^-2] and the outer cut-off radius [nm]
        o: output data dictionary
        j: harmonic fitted
        l: Fourier variable [nm]
        a: values of the Fourier amplitude
        m: function of the model fitted
        W: weigh the points

    Output:
        e: standard error of the fit (sigma^2)
    """
    am = m(*p, o, j, l) # model values
    if W:
        e = np.log(a)/l**2 - np.log(m(*p, o, j, l))/l**2
        i = np.arange(1, len(l)+1)
        w = len(l)*np.log((i+1)/i)/np.log(len(l)+1)
    else:
        e = a - am
        w = 1
    return np.sum(w*e**2)/(len(l)-2)

@beartype
def fits_data(
    m: ModelFunction,
    f: str,
    o: dict,
    **kwargs,
) -> dict:
    """
    Return information on the fits made of model m on output o.

    A fit is calculated for each harmonic and each interval of L.

    Input:
        m: model function
        f: lowpass filter name
        o: output data dictionary
      **d: initial density [nm^-2] (default: 5e14*1e-18)
      **r: initial outer cut-off radius [nm] (default: 200)
      **b: bounds for (d, r) (default: 5e10*1e-18->5e18*1e-18 / 1->np.inf)
      **j: restriction of harmonics (default: no restriction)

    Output:
        fitdat: fits data dictionary

    The quantities contained in fitdat are the following:
        'j' (ScalarList): harmonics
        'L' (ScalarList): maximum values ​​of L [nm]
        'd' (ScalarList): optimal density for the fit [nm^-2]
        'r' (ScalarList): optimal outer cut-off radius for the fit [nm]
        'e' (ScalarList): optimal fit errors
        'm' (Callable): model function
    """
    # optional parameters
    d = kwargs.pop('d', 5e14*1e-18) # initial density [nm^-2]
    r = kwargs.pop('r', 200) # initial outer cut-off radius [nm]
    b = kwargs.pop('b', ((5e10*1e-18, 5e18*1e-18), (1, np.inf))) # bounds
    j = kwargs.pop('j', o['j']) # harmonics to fit
    if len(kwargs)>0:
        raise ValueError("wrong keywords: "+str(kwargs))
    # fit
    J, L, D, R, E, = [], [], [], [], [] # fits information
    for h in j: # harmonic index
        i_j = o['index'][h] # harmonic index in outdat
        i_L = 3 # index of the maximum value of L for the fit
        failed = False # a fit has failed
        while not failed and i_L<=o[f][i_j]:
            fj = o['j'][i_j] # harmonic
            fl = o['L'][:i_L] # maximum value ​​of L
            fa = o['A'][i_j][:i_L] # Fourier amplitudes for harmonic j
            fr = scipy.optimize.minimize(
                standard_error,
                np.array((d, r)),
                args=(o, int(fj), fl, fa, m, f=='f2'),
                method='Nelder-Mead',
                bounds=b,
            )
            if fr.success:
                J.append(fj) # store harmonic
                L.append(fl[-1]) # store maximum value ​​of L [nm]
                D.append(fr.x[0]) # store optimal density [nm^-2]
                R.append(fr.x[1]) # store optimal outer cut-off radius [nm]
                E.append(fr.fun) # store error
            else:
                print(("/!\ "+fr.message+"\n    "
                    + "L="+str(fl[-1])+"nm "
                    + "j="+str(fj)+" "
                    + "model="+m.__name__
                ))
                failed = True
            i_L += 1
    fitdat = {
        'j': np.array(J), # harmonics
        'L': np.array(L), # maximum values of L [nm]
        'd': np.array(D), # optimal densities [nm^-2]
        'r': np.array(R), # optimal outer cut-off radii [nm]
        'e': np.array(E), # errors
        'm': m, # model function
    }
    return fitdat

@beartype
def plot(
    outdat: dict,
    **kwargs,
) -> None:
    """
    Export a figure with two representations of A(L) and the fits.

    Input:
        outdat: output data dictionary
      **expstm: export stem (default: output stem)
      **expdir: export directory (default: '')
      **exfmt: export format (default: 'pdf')
      **title: export title (default: output stem)
      **fitdat: fits data dictionary
      **j: restriction of harmonics (default: no restriction)
      **L: restriction of fits maximum L value (default: no restriction)
    """
    # optional parameters
    expstm = kwargs.pop('expstm', outdat['stm']) # export file stem
    expdir = kwargs.pop('expdir', '') # export directory
    expfmt = kwargs.pop('expfmt', 'pdf') # export format
    title = kwargs.pop('title', expstm.replace("_", " ")) # title
    fitdat = kwargs.pop('fitdat', None) # fits data dictionary
    j = kwargs.pop('j', outdat['j']) # harmonics restriction
    L = kwargs.pop('L', np.inf) # restriction of maximum L value
    if len(kwargs)>0:
        raise ValueError("wrong keywords: "+str(kwargs))
    # collect data
    data = [] # list of the information on the curves to display
    for h in j: # collect the output amplitudes to display
        i_j = outdat['index'][h] # index of the harmonic in c
        i_L = min(outdat['f1'][i_j]+5, len(outdat['L'])) # range to display
        data.append({
            'L': outdat['L'][:i_L], # x variable
            'A': outdat['A'][i_j][:i_L], # y variable
            'm': '.-' if fitdat is None else '.', # marker
            'n': "A_{"+str(h)+"}(L)", # label
            'c': 'C'+str(h-1), # color
            'l': '', # extra information
            })
    if not fitdat is None: # collect the fit amplitudes to display
        mask = np.in1d(fitdat['j'], j) # restrict to the requested harmonics
        if not L is None: # restrict the Fourier variable range
            mask = mask & np.in1d(fitdat['L'], L) # # restrict L
        maskj = fitdat['j'][mask] # apply mask on j
        maskL = fitdat['L'][mask] # apply mask on L
        maskd = fitdat['d'][mask] # apply mask on d
        maskr = fitdat['r'][mask] # apply mask on r
        maske = fitdat['e'][mask] # apply mask on e
        Lmin = outdat['L'][0] # minimum value of L (common to all fits)
        m = fitdat['m'] # model function
        for i in range(len(maskj)): # through remaining fits
            Lmax = maskL[i] # maximum value of L [nm]
            d, r = maskd[i], maskr[i] # density and outer cut-off radius
            v = r" \times 10^{".join(format(maske[i], '1.1e').split('e'))+"}"
            le = r"$ \hat{\sigma} ="+v+" $" # fit error
            ll = r"$ L \leq "+format(Lmax, '1.1f')+r" $" # fit range
            h = maskj[i] # harmonic
            l = np.linspace(Lmin, Lmax, 40) # range to plot the fit
            data.append({
                'L': l, # x variable
                'A': m(d, r, outdat, int(h), l), # y variable
                'm': '-', # marker
                'n': "A^{"+fitdat['m'].__name__+"}_{"+str(h)+r"}(L)", # label
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
    plt.savefig(os.path.join(expdir, expstm+'.'+expfmt), format=expfmt)
    plt.close('all')

@beartype
def export_model(
    modfun: ModelFunction,
    modflt: str,
    outdat: dict,
    **kwargs,
) -> None:
    """
    Export the information on fits and the corresponding figures.

    Input:
        modfun: model function
        modflt: model filter
        outdat: output data dictionary
      **expdir: export directory (default: '')
      **expstm: export stem (default: output stem)
      **fmtdat: data export format (default: 'dat')
      **fmtfit: fits export format (default: 'png')
      **title: title of the plots (default: model and output stem)
      **d: initial density [nm^-2] (default: see fits_data)
      **r: initial outer cut-off radius [nm] (default: see fits_data)
      **j: restriction of harmonics (default: see fits_data)
    """
    modnam = modfun.__name__ # model name
    # optional parameters
    expdir = kwargs.pop('expdir', '') # export directory
    expstm = kwargs.pop('expstm', outdat['stm']) # export file stem
    fmtdat = kwargs.pop('fmtdat', 'dat') # fits export format
    fmtfit = kwargs.pop('fmtfit', 'png') # fits export format
    title = kwargs.pop('title', modnam+" "+outdat['stm']) # fits title
    keys = [key for key in kwargs]
    kw = {key: kwargs.pop(key) for key in keys if key in ('d', 'r', 'j')}
    if len(kwargs)>0:
        raise ValueError("wrong keywords: "+str(kwargs))
    # collect
    fitdat = fits_data(modfun, modflt, outdat, **kw) # get fits data
    # export a figure for each fit
    moddir = expdir+"/fits_plot_"+modnam+"/" # model export directory
    if not os.path.exists(moddir):
        os.mkdir(moddir)
    for i in range(len(fitdat['e'])): # through the fits
        fitstm = ("j"
            + format(fitdat['j'][i], '01.0f')
            + "_"
            + format(fitdat['L'][i], '03.0f')
            + "nm"
        )
        plot(
            outdat, # output data
            expstm=fitstm, # stem of the fit
            expdir=moddir, # model export directory
            expfmt=fmtfit, # plot export format
            fitdat=fitdat, # fits
            title=title, # plot title
            j=np.array([fitdat['j'][i]]), # restriction of the harmonic
            L=np.array([fitdat['L'][i]]), # restriction of the maximum L
        )
    # export fits data
    fields = (
        'harmonic',
        'Lmax[nm]',
        'rho[m-2]',
        'Re[nm]',
        'error',
    )
    values = (
        fitdat['j'],
        fitdat['L'],
        fitdat['d']*1e18,
        fitdat['r'],
        fitdat['e'],
    )
    sep = " "
    fmt = ['%1.0f', '%5.1f', '%1e', '%1e', '%1e']
    with open(os.path.join(expdir, 'fits_data_'+modnam+"."+fmtdat), "w") as f:
        f.write(sep.join(fields)+'\n')
        np.savetxt(f, np.transpose(values), fmt=sep.join(fmt))

@beartype
def export(
    impstm: str,
    **kwargs,
) -> None:
    """
    Perform an analysis for each of the given models.

    Input:
        impstm: stem the simulation output (file or directory)
      **impdir: input directory (default: '')
      **expdir: export directory (default: '')
      **expstm: export stem (default: output stem)
      **fmtout: data export format (default: 'dat')
      **fmtdat: data export format (default: 'dat')
      **fmtfit: fits export format (default: 'png')
      **funflt: models and filters to analyze (default: GUW, WS, WC)
      **title: title of the plots (default: model and output stem)
      **d: initial density [nm^-2] (default: see export_model)
      **r: initial outer cut-off radius [nm] (default: see export_model)
      **j: restriction of harmonics (default: see export_model)
    """
    # optional parameters
    impdir = kwargs.pop('impdir', '') # import directory
    expdir = kwargs.pop('expdir', '') # export directory
    expstm = kwargs.pop('expstm', impstm+'_analysis') # export file stem
    fmtout = kwargs.pop('fmtout', 'pdf') # fits export format
    funflt = kwargs.pop('funflt', [
        (models.GUW2, 'f2'),
        (models.W1, 'f1'),
        (models.W2, 'f2'),
    ])
    keys = [key for key in kwargs]
    kw1 = {
        key: kwargs.pop(key) for key in keys
        if key in ('d', 'r', 'j', 'fmtfit', 'fmtdat', 'title')
    }
    kw2 = {key: val for key, val in kw1.items() if key in ('title', 'j')}
    if len(kwargs)>0:
        raise ValueError("wrong keywords: "+str(kwargs))
    # export directory
    dirstm = os.path.join(expdir, expstm)
    if not os.path.exists(dirstm):
        os.mkdir(dirstm)
    # load output data
    outdat = output_data(impstm, impdir=impdir)
    # plot output data
    plot(
        outdat,
        expstm='output_plot',
        expdir=dirstm,
        expfmt=fmtout,
        **kw2,
    )
    # fits
    for modfun, modflt in funflt:
        export_model(
            modfun,
            modflt,
            outdat,
            expdir=dirstm,
            expstm=expstm,
            **kw1
        )
