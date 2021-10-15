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

defmodspe = (
    (
        models.GUW1,
        'f1',
        np.array((5e14*1e-18, 1000, 5e-4, 10)),
        ((1e13*1e-18, 1e16*1e-18), (1, np.inf), (0, np.inf), (1, np.inf)),
    ),
    (
        models.GUW2,
        'f2',
        np.array((5e14*1e-18, 1000)),
        ((1e13*1e-18, 1e16*1e-18), (1, np.inf)),
    ),
    (
        models.W1,
        'f1',
        np.array((5e14*1e-18, 1000)),
        ((1e13*1e-18, 1e16*1e-18), (1, np.inf)),
    ),
    (
        models.W2,
        'f2',
        np.array((5e14*1e-18, 1000)),
        ((1e13*1e-18, 1e16*1e-18), (1, np.inf)),
    ),
)

matplotlib.use("Agg") # to export plots with no bitmap allocation errors

@beartype
def output_data(
    impstm: str,
    **kwargs,
) -> dict:
    """
    Return a set of quantities related to a simulation output.

    Input:
        impstm (str): stem of the simulation output
      **impdir (str): import directory (default: see collect.load_file)
      **frrprt (Callable): Fourier transform part (default: np.real)

    Output:
        outdat (dict): output data dictionary

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
        'frrprt' (Callable): Fourier transform part
    """
    # optional parameters
    frrprt = getkwa('frrprt', kwargs, Callable, np.real)
    # collect data
    outdat = {'stm': impstm} # output data dictionary
    qtynam = ['g', 'z', 'b', 'C', 'a', 'L', 'A'] # quantities to load
    for key, val in zip(qtynam, collect.load(qtynam, impstm, **kwargs)):
        outdat[key] = val # store the loaded quantities
    outdat['A'] = frrprt(outdat['A'])
    outdat['g'] = outdat['g']/outdat['a'] # correct diffraction vector norm
    outdat['j'] = np.arange(len(outdat['A'])) + 1 # harmonics
    # pre-calculated quantities
    vertij = outdat['j'].reshape((len(outdat['j']), 1))
    outdat['jg'] = vertij*outdat['g'] # j*g [nm^-1]
    outdat['jgb'] = np.sum(outdat['b']*outdat['jg'], axis=1) # (j*g).b [1]
    outdat['b2'] = np.sum(outdat['b']**2) # |b|^2 [nm^2]
    outdat['jg2'] = np.sum(outdat['jg']**2, axis=1) # |j*g|^2 [nm^-2]
    # filters
    outdat['f0'] = [filters.f0(a) for a in outdat['A']] # suppress noise
    outdat['f1'] = [filters.f1(a) for a in outdat['A']] # suppress noise
    outdat['f2'] = [filters.f2(a, outdat['L']) for a in outdat['A']] # linear
    # reverse index directory
    outdat['index'] = {outdat['j'][i]: i for i in range(len(outdat['j']))}
    outdat['frrprt'] = frrprt
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
        p (ScalarList): optimization parameters
        o (dict): output data dictionary
        j (int): harmonic fitted
        l (ScalarList): Fourier variable [nm]
        a (ScalarList): values of the Fourier transform
        m (ModelFunction): function of the model fitted
        W (bool): weigh the points

    Output:
        e (Scalar): standard error of the fit (sigma^2)
    """
    am = m(p, o, j, l) # model values
    if W:
        e = np.log(a)/l**2 - np.log(m(p, o, j, l))/l**2
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
    p: ScalarList,
    **kwargs,
) -> dict:
    """
    Return information on the fits made of model m on output o.

    A fit is calculated for each harmonic and each interval of L.

    Input:
        m (ModelFunction): model function
        f (str): lowpass filter name
        o (dict): output data dictionary
        p (ScalarList): initial optimization parameters
      **b (tuple): initial optimization parameters bounds (default: None)
      **j (ScalarList): restriction of harmonics (default: no restriction)

    Output:
        fitdat (dict): fits data dictionary

    The quantities contained in fitdat are the following:
        'j' (ScalarList): harmonics
        'l' (ScalarList): maximum values ​​of L [nm]
        'e' (ScalarList): optimal fit standard errors
        'p' (ScalarList): optimal parameters
        'm' (Callable): model function
    """
    # optional parameters
    b = getkwa('b', kwargs, Optional[tuple], None)
    j = getkwa('j', kwargs, ScalarList, o['j'])

    endkwa(kwargs)
    # fit
    J, L, E, P = [], [], [], [] # fits information
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
                p,
                args=(o, int(fj), fl, fa, m, f=='f2'),
                method='Nelder-Mead',
                bounds=b,
            )
            if fr.success:
                J.append(fj) # store harmonic
                L.append(fl[-1]) # store maximum value ​​of L [nm]
                E.append(fr.fun) # store error
                p = fr.x
                P.append(p) # store optimal parameters [nm^-2]
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
        'l': np.array(L), # maximum values of L [nm]
        'e': np.array(E), # errors
        'p': np.array(P), # optimal parameters
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
        outdat (str): output data dictionary
      **expstm (str): export stem (default: output stem)
      **expdir (str): export directory (default: '')
      **expfmt (str): export format (default: 'pdf')
      **figttl (str): figure title (default: based on expstm)
      **fitdat (NoneType|dict): fits data dictionary
      **j (ScalarList): restriction of harmonics (default: no restriction)
      **L (Scalar): restriction of fits maximum L value (default: np.inf)
    """
    # optional parameters
    expstm = getkwa('expstm', kwargs, str, outdat['stm'])
    expdir = getkwa('expdir', kwargs, str, '')
    expfmt = getkwa('expfmt', kwargs, str, 'pdf')
    figttl = getkwa('figttl', kwargs, str, expstm.replace("_", " "))
    fitdat = getkwa('fitdat', kwargs, Optional[dict], None)
    j = getkwa('j', kwargs, ScalarList, outdat['j'])
    L = getkwa('L', kwargs, Scalar, np.inf)
    endkwa(kwargs)
    # Fourier transform part
    if outdat['frrprt'] == np.real:
        prtnot = r"\mathcal{R}(FUN)" # notation for real part
    elif outdat['frrprt'] == np.absolute:
        prtnot = r"|FUN|" # notation for module
    else:
        prtnot = r"\mathrm{"+frrprt.__name__+r"}(FUN)"
    # collect data
    data = [] # list of the information on the curves to display
    for h in j: # collect the output amplitudes to display
        i_j = outdat['index'][h] # index of the harmonic in c
        i_L = min(outdat['f1'][i_j], len(outdat['L'])) # range to display
        data.append({
            'L': outdat['L'][:i_L], # x variable
            'A': outdat['A'][i_j][:i_L], # y variable
            'm': '.-' if fitdat is None else '.', # marker
            'n': prtnot.replace("FUN", "A_{"+str(h)+"}(L)"), # label
            'c': 'C'+str(h-1), # color
            'l': '', # extra information
            })
    if not fitdat is None: # collect the fit amplitudes to display
        mask = np.in1d(fitdat['j'], j) # restrict to the requested harmonics
        if not L is None: # restrict the Fourier variable range
            mask = mask & np.in1d(fitdat['l'], L) # # restrict L
        maskj = fitdat['j'][mask] # apply mask on j
        maskl = fitdat['l'][mask] # apply mask on L
        maske = fitdat['e'][mask] # apply mask on e
        maskp = fitdat['p'][mask] # apply mask on d
        Lmin = outdat['L'][0] # minimum value of L (common to all fits)
        m = fitdat['m'] # model function
        for i in range(len(maskj)): # through remaining fits
            Lmax = maskl[i] # maximum value of L [nm]
            p = maskp[i] # density and outer cut-off radius
            v = r" \times 10^{".join(format(maske[i], '1.1e').split('e'))+"}"
            le = r"$ \hat{\sigma} ="+v+" $" # fit error
            ll = r"$ L \leq "+format(Lmax, '1.1f')+r" $" # fit range
            h = maskj[i] # harmonic
            l = np.linspace(Lmin, Lmax, 40) # range to plot the fit
            data.append({
                'L': l, # x variable
                'A': m(p, outdat, int(h), l), # y variable
                'm': '-', # marker
                'n': fitdat['m'].__name__+"_{"+str(h)+r"}(L)", # label
                'c': 'C'+str(h-1), # color
                'l': ", "+le+", "+ll, # extra information
            })
    # fig
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.1)
    fig.suptitle(figttl)
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
            label=r"$\ln("+d['n']+r")/L^2$"+d['l'],
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
    p:  ScalarList,
    **kwargs,
) -> None:
    """
    Export the information on fits and the corresponding figures.

    Input:
        modfun (Callable): model function
        modflt (str): model filter name
        outdat (dict): output data dictionary
        p (ScalarList): initial optimization parameters
      **expdir (str): export directory (default: '')
      **expstm (str): export stem (default: output stem)
      **fmtdat (str): data export format (default: 'dat')
      **fmtfit (str): fits export format (default: 'png')
      **figttl (str): title of the plots (default: model and output stem)
      **valsep (str): value separator (default: ' ')
      **prmfmt (tuple): optimization parameter format (default: ('1e', ...))
      **prmnam (tuple): optimization parameter names (default: ('p1', ...))
      **j (ScalarList): restriction of harmonics (default: see fits_data)
      **b (tuple): initial optimization parameters bounds (default: None)
    """
    modnam = modfun.__name__ # model name
    # optional parameters
    expdir = getkwa('expdir', kwargs, str, '')
    expstm = getkwa('expstm', kwargs, str, outdat['stm'])
    fmtdat = getkwa('fmtdat', kwargs, str, 'dat')
    fmtfit = getkwa('fmtfit', kwargs, str, 'png')
    figttl = getkwa('figttl', kwargs, str, modnam+" "+outdat['stm'])
    valsep = getkwa('valsep', kwargs, str, ' ')
    prmfmt = getkwa('prmfmt', kwargs, tuple, tuple(['%1e' for i in p]))
    prmnam = getkwa('prmnam', kwargs, tuple,
        tuple(['p'+str(i+1) for i in range(len(p))]))
    # collect
    fitdat = fits_data(modfun, modflt, outdat, p, **kwargs) # get fits data
    # export a figure for each fit
    moddir = expdir+"/fits_plot_"+modnam+"/" # model export directory
    if not os.path.exists(moddir):
        os.mkdir(moddir)
    for i in range(len(fitdat['e'])): # through the fits
        fitstm = ("j"+ format(fitdat['j'][i], '01.0f')
                + "_"+ format(fitdat['l'][i], '03.0f')+ "nm")
        plot(
            outdat, # output data
            expstm=fitstm, # stem of the fit
            expdir=moddir, # model export directory
            expfmt=fmtfit, # plot export format
            fitdat=fitdat, # fits
            figttl=figttl, # figure title
            j=np.array([fitdat['j'][i]]), # restriction of the harmonic
            L=fitdat['l'][i], # restriction of the maximum L
        )
    # export fits data
    fields = ['j', 'Lmax', 'error']+list(prmnam)
    fmt = ['%1.0f', '%5.1f', '%1e']+list(prmfmt)
    values = np.concatenate((
        np.transpose((fitdat['j'], fitdat['l'], fitdat['e'])),
        fitdat['p'],
    ), axis=1)
    with open(os.path.join(expdir, 'fits_data_'+modnam+"."+fmtdat), "w") as f:
        f.write(valsep.join(fields)+'\n')
        np.savetxt(f, values, fmt=valsep.join(fmt))

@beartype
def export(
    impstm: str,
    **kwargs,
) -> None:
    """
    Perform an analysis for each of the given models.

    Input:
        impstm (str): stem the simulation output (file or directory)
      **impdir (str): input directory (default: '')
      **expdir (str): export directory (default: '')
      **expstm (str): export stem (default: output stem)
      **fmtout (str): data export format (default: 'dat')
      **fmtdat (str): data export format (default: 'dat')
      **fmtfit (str): fits export format (default: 'png')
      **modspe (tuple): models specifications (default: defmodspe)
      **figttl (str): title of the plots (default: model and output stem)
      **j (str): restriction of harmonics (default: see export_model)
    """
    # optional parameters
    impdir = getkwa('impdir', kwargs, str, '')
    expdir = getkwa('expdir', kwargs, str, '')
    expstm = getkwa('expstm', kwargs, str, impstm+'_analysis')
    fmtout = getkwa('fmtout', kwargs, str, 'pdf')
    modspe = getkwa('modspe', kwargs, tuple, defmodspe)
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
        **{key: val for key, val in kwargs.items() if key in ('figttl', 'j')},
    )
    # fits
    for modfun, modflt, p, b in modspe:
        export_model(
            modfun,
            modflt,
            outdat,
            p,
            b=b,
            expdir=dirstm,
            expstm=expstm,
            **kwargs
        )
