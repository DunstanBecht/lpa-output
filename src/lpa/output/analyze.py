#!/usr/bin/env python
# coding: utf-8

"""
Tools for the analysis of the simulation output.
"""

import matplotlib.pyplot as plt
import scipy.optimize
from . import *
from . import collect
from . import models

@beartype
def A(
    n: str,
    j: ScalarList,
    p: str = 'output/',
) -> ScalarListList:
    """
    Return the averaged Fourier amplitude for the harmonics j from n.

    Input:
        n: name of the results directory of the analyzed sample
        j: harmonics of diffraction vector
        p: path where n can be found

    Output:
        a: averaged Fourier amplitude for each harmonic in j
        [
            np.array([A_1(L_1), A_1(L_2), A_1(L_3), ...]),
            np.array([A_2(L_1), A_2(L_2), A_2(L_3), ...]),
            ...
        ]
    """
    q = []
    for h in j:
        if h == 1:
            q += ['cos_AL', 'sin_AL']
        else:
            q += ['Cos_'+str(h)+'AL', 'Sin_'+str(h)+'AL']
    m = collect.load(q, n, p)
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
    s: str,
    p: str = 'output/',
    j: Optional[ScalarList] = None,
) -> dict:
    """
    Return a set of common quantities related to a simulation.

    Input:
        s: name of the results directory of the analyzed sample
        p: path where s can be found
        j: harmonics studied

    Output:
        c: dictionary containing the quantities common to a simulation

    The quantities contained in c are the following:
        'name': name of the results directory
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
        'index' (dict): index of the harmonics
    """
    c = {}
    q = ['g', 'z', 'b', 'C', 'a', 'J', 'L']
    for k, v in zip(q, collect.load(q, s, p)):
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
    c['A'] = A(s, c['j'], p)
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
) -> dict:
    """
    Return the table of fits of model m on simulation c.

    A fit is calculated for each harmonic and each interval of L.

    Input:
        m: model function
        c: dictionary containing the quantities common to a simulation
        f: lowpass filter name

    Output:
        n: column names
        v: column values

    The table contains the following columns:
        j: harmonics
        L: maximum values ​​of L [nm]
        d: optimal density for the fit [nm^-2]
        r: optimal outer cut-off radius for the fit [nm]
        e: fit errors
    """
    J, L, D, R, E, = [], [], [], [], []
    for i_j in range(len(c['j'])):
        for i_L in range(3, c[f][i_j]):
            j = c['j'][i_j]
            l = c['L'][:i_L]
            a = c['A'][i_j][:i_L]
            def error(p)-> Scalar:
                return np.sum((a-m(*p, c, int(j), l))**2/(i_L-2))
            print('')
            p = scipy.optimize.fmin(error, (0.02, 200), ftol=1e-10)
            print('')
            J.append(j)
            L.append(l[-1])
            D.append(p[0])
            R.append(p[1])
            E.append(error(p))
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
    n: str,
    p: str,
    j: Optional[ScalarList] = None,
    t: str = "",
    f: dict = None,
    L: Optional[ScalarList] = None,
    e: str = "pdf",
) -> None:
    """
    Export a figure with two representations of A(L) and the fits.

    Input:
        c: common parameters dictionary
        n: name of the exported figure
        p: path where the figure must be exported
        j: restriction of harmonics to be displayed
        t: title of the figure
        f: models fits
        L: restriction of maximum L values
        e: file extension
    """
    if j is None:
        j = c['j']
    data = []
    for h in j:
        i_j = c['index'][h]
        i_L = c['i1'][i_j]
        data.append({
            'L': c['L'][:i_L],
            'A': c['A'][i_j][:i_L],
            'm': '.',
            'n': "A_{"+str(h)+"}(L)",
            'c': 'C'+str(h-1),
            'l': '',
            })

    if f != None:
        mask = np.in1d(f['j'], j)
        if L != None:
            mask = mask & np.in1d(f['L'], L)

        maskj = f['j'][mask]
        maskL = f['L'][mask]
        maskd = f['d'][mask]
        maskr = f['r'][mask]
        maske = f['e'][mask]
        min = c['L'][0]
        m = f['m']
        for i in range(len(maskj)):
            max = maskL[i]
            d, r = maskd[i], maskr[i]
            v = r" \times 10^{".join(format(maske[i], '1.1e').split('e'))+"}"
            le = r"$ \sigma^2 ="+v+" $"
            ll = r"$ L \leq "+str(max)+r" $"
            h = maskj[i]
            l = np.linspace(min, max, 20)
            data.append({
                'L': l,
                'A': m(d, r, c, int(h), l),
                'm': '-',
                'n': "A^"+f['name'][0]+"_{"+str(h)+r"}(L)",
                'c': 'C'+str(h-1),
                'l': ", "+le+", "+ll,
            })
    # fig
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.1)
    fig.suptitle(t)
    # ax1
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
    # ax2
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
    plt.savefig(p+n+'.'+e)
    plt.close('all')

@beartype
def export_model(
    m: dict,
    c: dict,
    p: str,
    n: str,
) -> None:
    """
    Export the table of fits and the corresponding figures.

    Input:
        m: dictionary containing information about the model
        c: common quantities
        p: path
        n: file name
    """
    f = fit(m['function'], c, m['filter'])
    p_f = p+"/"+f['name']+"/"
    if not os.path.exists(p_f):
        os.mkdir(p_f)
    for i in range(len(f['j'])):
        pn = ("j"
            + format(f['j'][i], '01.0f')
            + "_"
            + format(f['L'][i], '03.0f')
            + "nm"
        )
        plot(
            c,
            pn,
            p_f,
            f=f,
            j=np.array([f['j'][i]]),
            L=np.array([f['L'][i]]),
            e='png',
        )
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
    with open(p+n+"_"+f['name']+".csv", "w") as f:
        f.write(sep.join(fields)+'\n')
        np.savetxt(f, np.transpose(values), fmt=sep.join(fmt))

@beartype
def export(
    s: str,
    i: str = "",
    o: str = "",
) -> None:
    """
    Do an analysis of A with the available models.

    Input:
        s: name of the results directory of the analyzed sample
        i: path where s can be found
        o: path where to export the analysis
    """
    if i!="" and i[-1]!="/":
        i += "/"
    if o!="" and o[-1]!="/":
        o += "/"
    # title
    t = " ".join(s.split("_"))
    # directory
    o_s = o+s+"/"
    if not os.path.exists(o_s):
        os.mkdir(o_s)
    # general
    c = common(s, i, np.array([1, 2]))
    plot(c, s, o_s, t=t)
    # models
    analyzed = [
        {'function': models.Groma, 'filter': 'i2',},
        {'function': models.Kamminga, 'filter': 'i2',},
        {'function': models.Wilkens, 'filter': 'i1',},
    ]
    for m in analyzed:
        export_model(m, c, o_s, s)
