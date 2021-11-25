#!/usr/bin/env python
# coding: utf-8

"""
Tools for the analysis of the simulation output.
"""

import matplotlib
import matplotlib.pyplot as plt
import scipy.optimize
from . import *
from . import __version__
from . import collect, models, filters

matplotlib.use("Agg") # to export plots with no bitmap allocation errors

@beartype
def defmodspe(
    d: Scalar,
    r: Scalar = 1000,
) -> tuple:
    """
    Return a default model specifications.

    Input:
        d (Scalar): dislocation density initial guess [nm^-2]
      **r (Scalar): outer cut-off radius initial guess [nm^-2]

    Output:
        modspe (tuple): model specifications
    """
    modspe = (
        (
            models.GUW1,
            'F1',
            np.array((d, r, 0, 600)),
            ('density[nm-2]', 'cut-off-radius[nm]', 'fluctuation', 'R0[nm]'),
            ((1e-18, 1), (1e-20, 1e5), (0, 1e2), (1e-200, 1e5)),
        ),
        (
            models.GUW2,
            'F2',
            np.array((d, r)),
            ('density[nm-2]', 'cut-off-radius[nm]'),
            ((1e-18, 1), (1e-20, 1e5)),
        ),
        (
            models.W1,
            'F1',
            np.array((d, r)),
            ('density[nm-2]', 'cut-off-radius[nm]'),
            ((1e-18, 1), (1e-20, 1e5)),
        ),
        (
            models.W2,
            'F2',
            np.array((d, r)),
            ('density[nm-2]', 'cut-off-radius[nm]'),
            ((1e-18, 1), (1e-20, 1e5)),
        ),
    )
    return modspe

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
        'd' (Scalar): dislocation density [nm^-2]
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
        'F1' (List): index at which the noise starts for each harmonic
        'F2' (List): index of linear zone end for each harmonic
        'index' (dict): reverse indexing for the harmonics
        'frrprt' (Callable): Fourier transform part
    """
    # optional parameters
    frrprt = getkwa('frrprt', kwargs, Callable, np.real)
    # collect data
    outdat = {'stm': impstm} # output data dictionary
    qtynam = ['A', 'd', 'g', 'z', 'b', 'C', 'a', 'L'] # quantities to load
    for key, val in zip(qtynam, collect.load(qtynam, impstm, **kwargs)):
        outdat[key] = val # store the loaded quantities
    outdat['A'] = frrprt(outdat['A'])
    outdat['d'] = outdat['d']*1e-18 # change density unit
    outdat['g'] = outdat['g']/outdat['a'] # correct diffraction vector norm
    outdat['b'] = outdat['b']*outdat['a']/2 # correct Burgers vector norm
    outdat['j'] = np.arange(len(outdat['A'])) + 1 # harmonics
    # pre-calculated quantities
    vertij = outdat['j'].reshape((len(outdat['j']), 1))
    outdat['jg'] = vertij*outdat['g'] # j*g [nm^-1]
    outdat['jgb'] = np.sum(outdat['b']*outdat['jg'], axis=1) # (j*g).b [1]
    outdat['b2'] = np.sum(outdat['b']**2) # |b|^2 [nm^2]
    outdat['jg2'] = np.sum(outdat['jg']**2, axis=1) # |j*g|^2 [nm^-2]
    # filters
    outdat['F0'] = [filters.F0(a) for a in outdat['A']] # suppress noise
    outdat['F1'] = [filters.F1(a) for a in outdat['A']] # suppress noise
    outdat['F2'] = [filters.F2(a, outdat['L']) for a in outdat['A']] # linear
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
        e (Scalar): standard error of the fit
    """
    am = m(p, o, j, l) # model values
    if W:
        e = np.log(a/am)/l**2
        i = np.arange(1, len(l)+1)
        w = len(l)*np.log((i+1)/i)/np.log(len(l)+1)
    else:
        e = a - am
        w = 1
    return np.sqrt(np.sum(w*e**2)/(len(l)-2))

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
      **b (tuple): initial optimization parameter bounds (default: defmodspe)
      **j (ScalarList): restriction of harmonics (default: no restriction)
      **maxitr (int): maximum number of iterations (default: 1e7)

    Output:
        fitdat (dict): fits data dictionary

    The quantities contained in fitdat are the following:
        'j' (ScalarList): harmonics
        'l' (ScalarList): maximum values ​​of L [nm]
        'e' (ScalarList): optimal fit standard errors
        'p' (ScalarList): optimal parameters
        'm' (Callable): model function
        'f' (str): filter used
    """
    # optional parameters
    b = getkwa('b', kwargs, Optional[tuple], None)
    j = getkwa('j', kwargs, ScalarList, o['j'])
    maxitr = getkwa('maxitr', kwargs, int, 10**7)
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
            fr = lambda prm: scipy.optimize.minimize(
                standard_error,
                prm,
                args=(o, int(fj), fl, fa, m, f=='F2'),
                method='Nelder-Mead',
                bounds=b,
                options={'maxiter': maxitr}
            )
            tstprm = [p] # tested initial guess parameters
            if len(P)>0:
                tstprm.append(P[-1])
            tstres = [fr(prm) for prm in tstprm]
            tstres = [res for res in tstres if res.success]
            if len(tstres)>0:
                tstres = sorted(tstres, key=lambda res: res.fun)
                res = tstres[0]
                J.append(fj) # store harmonic
                L.append(fl[-1]) # store maximum value ​​of L [nm]
                E.append(res.fun) # store error
                P.append(res.x) # store optimal parameters [nm^-2]
                print((f"dst: {res.x[0]/o['d']:11.9f} | "
                       f"rad: {res.x[1]:9.3e} | "
                       f"err: {res.fun:9.3e} | "
                       f"{m.__name__:>4}_{fj} | "
                       f"p3+: {res.x[2:]}"))
            else:
                print((f"/!\ {' & '.join([res.message for res in tstres])}\n "
                       f"L={fl[-1]}nm "
                       f"j={fj} "
                       f"model={m.__name__}"))
                failed = True
            i_L += 1
    fitdat = {
        'j': np.array(J), # harmonics
        'l': np.array(L), # maximum values of L [nm]
        'e': np.array(E), # errors
        'p': np.array(P), # optimal parameters
        'm': m, # model function
        'f': f, # filter used
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
        prtnot = fr"\mathrm{{ {frrprt.__name__} }}(FUN)"
    # collect data
    data = [] # list of the information on the curves to display
    for h in j: # collect the output amplitudes to display
        i_j = outdat['index'][h] # index of the harmonic in c
        i_L = min(outdat['F1'][i_j], len(outdat['L'])) # range to display
        data.append({
            'L': outdat['L'][:i_L], # x variable
            'A': outdat['A'][i_j][:i_L], # y variable
            'm': '.-' if fitdat is None else '.', # marker
            'n': prtnot.replace("FUN", fr"A_{{ {h} }}(L)"), # label
            'c': f'C{h-1}', # color
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
            v = format(maske[i], '1.1e').replace('e', r"\times 10 ^{")+"}"
            le = fr"$ \hat{{\sigma}} = {v} $" # fit error
            ll = fr"$ L \leq {Lmax:1.1f} $" # fit range
            h = maskj[i] # harmonic
            l = np.linspace(Lmin, Lmax, 40) # range to plot the fit
            data.append({
                'L': l, # x variable
                'A': m(p, outdat, int(h), l), # y variable
                'm': '-', # marker
                'n': f"{fitdat['m'].__name__}_{{ {h} }}(L)", # label
                'c': f'C{h-1}', # color
                'l': f", {le}, {ll}", # extra information
            })
    # fig
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.subplots_adjust(left=0.06, right=0.95, bottom=0.1)
    fig.suptitle(figttl)
    # ax1: A(L) as a function of L (log scale)
    for d in data:
        ax1.plot(
            d['L'], # x variable
            d['A'], # y variable
            d['m'], # marker
            label=fr"$ {d['n']} $ {d['l']}",
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
            label=fr"$ \ln({d['n']}) / L^2 $ {d['l']}",
            color=d['c'],
        )
    ax2.set_xscale("log")
    ax2.set_xlabel(r"$L \ (nm)$")
    ax2.grid()
    ax2.legend()
    # version
    ax2.text(
        1.05,
        0.5,
        f"lpa-output ({__version__})",
        rotation=90,
        fontfamily='monospace',
        ha='left',
        va='center',
        transform=ax2.transAxes,
    )
    # fit representation
    if not fitdat is None:
        ax = {'F1': ax1, 'F2': ax2}[fitdat['f']]
        for spine in ax.spines.values():
            spine.set_edgecolor('C3')
        ax.set_title(
            f"(representation used to calculate the fitting deviation)",
            color='C3',
        )
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
      **prmnam (tuple): optimization parameter names (default: ('p1', ...))
      **maxitr (int): maximum number of iterations (default: 1e7)
      **j (ScalarList): restriction of harmonics (default: see fits_data)
      **b (tuple): initial optimization parameter bounds (default: None)
    """
    modnam = modfun.__name__ # model name
    # optional parameters
    expdir = getkwa('expdir', kwargs, str, '')
    expstm = getkwa('expstm', kwargs, str, outdat['stm'])
    fmtdat = getkwa('fmtdat', kwargs, str, 'dat')
    fmtfit = getkwa('fmtfit', kwargs, str, 'png')
    figttl = getkwa('figttl', kwargs, str, modnam+" "+outdat['stm'])
    valsep = getkwa('valsep', kwargs, str, ' ')
    prmnam = getkwa('prmnam', kwargs, tuple,
        tuple(['p'+str(i+1) for i in range(len(p))]))
    # collect
    fitdat = fits_data(modfun, modflt, outdat, p, **kwargs) # get fits data
    # export a figure for each fit
    moddir = f"{expdir}/fits_plot_{modnam}/" # model export directory
    if not os.path.exists(moddir):
        os.mkdir(moddir)
    for i in range(len(fitdat['e'])): # through the fits
        plot(
            outdat, # output data
            expstm=f"j{fitdat['j'][i]:01.0f}_{fitdat['l'][i]:03.0f}nm",
            expdir=moddir, # model export directory
            expfmt=fmtfit, # plot export format
            fitdat=fitdat, # fits
            figttl=figttl, # figure title
            j=np.array([fitdat['j'][i]]), # restriction of the harmonic
            L=fitdat['l'][i], # restriction of the maximum L
        )
    # export fits data
    dftcol = f"# j{valsep}Lmax[nm]{valsep}   error{valsep}"
    dftfmt = f"%3.0f{valsep}%8.1f{valsep}%8.1e{valsep}"
    col = dftcol + valsep.join([format(n, '>22') for n in prmnam])
    fmt = dftfmt + valsep.join(['%22.15e' for n in prmnam])
    val = np.concatenate((
        np.transpose((fitdat['j'], fitdat['l'], fitdat['e'])),
        fitdat['p'],
    ), axis=1)
    with open(os.path.join(expdir, f'fits_data_{modnam}.{fmtdat}'), "w") as f:
        f.write(f"{__version__:>8} # v: lpa-ouput version\n")
        f.write(f"{outdat['frrprt'].__name__:>8} # t: A(L) transformation\n")
        f.write(f"{modfun.__name__:>8} # m: model function\n")
        f.write(f"{modflt:>8} # f: filter\n")
        f.write(col+'\n')
        np.savetxt(f, val, fmt=fmt)

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
      **figttl (str): title of the plots (default: model and output stem)
      **frrprt (Callable): Fourier transform part (default: np.real)
      **modspe (tuple): models specifications (default: defmodspe)
      **maxitr (int): maximum number of iterations (default: 1e7)
      **modspe (tuple): model specifications (default: defmodspe)
      **j (str): restriction of harmonics (default: see export_model)
    """
    # optional parameters
    impdir = getkwa('impdir', kwargs, str, '')
    expdir = getkwa('expdir', kwargs, str, '')
    expstm = getkwa('expstm', kwargs, str, impstm+'_analysis')
    fmtout = getkwa('fmtout', kwargs, str, 'pdf')
    frrprt = getkwa('frrprt', kwargs, Callable, np.real)
    # export directory
    dirstm = os.path.join(expdir, expstm)
    if not os.path.exists(dirstm):
        os.mkdir(dirstm)
    # load output data
    outdat = output_data(impstm, impdir=impdir, frrprt=frrprt)
    modspe = getkwa('modspe', kwargs, tuple, defmodspe(outdat['d']))
    # plot output data
    plot(
        outdat,
        expstm='output_plot',
        expdir=dirstm,
        expfmt=fmtout,
        **{key: val for key, val in kwargs.items() if key in ('figttl', 'j')},
    )
    # fits
    for modfun, modflt, p, prmnam, b in modspe:
        export_model(
            modfun,
            modflt,
            outdat,
            p,
            prmnam=prmnam,
            b=b,
            expdir=dirstm,
            expstm=expstm,
            **kwargs
        )
