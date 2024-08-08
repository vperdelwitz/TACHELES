from PyAstronomy import pyasl
from astropy.io import fits
from math import pi
from numba import njit
from numpy.polynomial.polynomial import Polynomial
from photutils.aperture import CircularAperture, aperture_photometry
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.interpolate import interpn
import batman
import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import os.path
import scipy
import scipy.stats
import time

def lc_chrom_calc(t, Rp, logH, dlogH, logHup, logHlow, Rmax, Ns, period, a, inc, e, t0, w, Omega, npl, osf):
    chrom_r, logHa, chrom_tab = get_chrom_tab(Rmax,Ns,dlogH,logHup,logHlow)
    chrom_Fr = interp1d(logHa, chrom_tab, axis=0)(logH)
    Fr_fun = interp1d(chrom_r, chrom_Fr,fill_value='extrapolate')
    xyp = pyasl.KeplerEllipse(a, period, e=e, tau=t0, Omega=Omega, i=inc, w=w).xyzPos(t)
    xpm, ypm, apm = makepl(npl, osf, plotting=0)
    lc_chrom = calc_lc_from_rprof(xyp, Rp, Fr_fun, Ns, xpm, ypm, apm)
    return lc_chrom

def build_chrom_tab(logHa, Rmax, Ns):
    print("Calculating chrom table...")
    chrom_tab = np.zeros((logHa.size, Ns*Rmax))
    for iH, logH in enumerate(logHa):
        print('Working on: ',iH,'/',logHa.size,'logH=',logH,end="\r")
        chrom_r, chrom_tab[iH, :] = calc_chrom_rprof(10**logH, Rmax, Ns)
    return chrom_r, chrom_tab

def get_chrom_tab(Rmax, Ns, dlogH, logHup, logHlow, recalculate=False):
    logHa=np.arange(logHlow,logHup,dlogH)
    Htab_fn=(f'Htab_Rmax{Rmax}_Ns{Ns}_dlogH{dlogH}.mat')
    if not os.path.exists(Htab_fn) or recalculate:
        print('Preparing ', Htab_fn)
        chrom_r, chrom_tab = build_chrom_tab(logHa, Rmax, Ns)
        scipy.io.savemat(Htab_fn,{'logHa':logHa,'chrom_r':chrom_r,'chrom_tab':chrom_tab,'Rmax':Rmax,'Ns':Ns})
    else:
        print('Reading ', Htab_fn)
        dic=scipy.io.loadmat(Htab_fn)
        for key, value in dic.items():
            logHa=np.squeeze(dic['logHa'])
            chrom_r=np.squeeze(dic['chrom_r'])
            chrom_tab=np.squeeze(dic['chrom_tab'])

    return chrom_r, logHa, chrom_tab

def interpolate_chrom_rprof(logH, logHa, chrom_tab):
    chrom_Fr = interp1d(logHa, chrom_tab, axis=0)(logH)
    return chrom_Fr

def rebin(a, shape):
    sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def makepl(npl, osf, plotting=0):  # npl = number of planet pixels, osf = oversample factor
    xpl, ypl = np.meshgrid(np.linspace(-1.1, 1.1, npl * osf), np.linspace(-1.1, 1.1, npl * osf))
    rpl = np.sqrt(np.square(xpl) + np.square(ypl))
    apl = np.ones((npl * osf, npl * osf), dtype='double')
    apl[rpl > 1] = 0
    if osf>1:
        apl = rebin(apl, (npl, npl))
        xpl = rebin(xpl, (npl, npl))
        ypl = rebin(ypl, (npl, npl))
        rpl = rebin(rpl, (npl, npl))

    if plotting:
        plt.imshow(apl)
        plt.title('Planet')
        plt.colorbar()
        plt.show()

    wpl = np.nonzero(rpl <= 1)
    xpl = xpl[wpl]
    ypl = ypl[wpl]
    apl = apl[wpl]
    apl = apl * pi / apl.sum()  # area of a planet with radius 1 is pi

    return xpl, ypl, apl

def calc_chrom_rprof(H, Rmax, Ns, plotting=0):
    def integrand(y, x, H):
        r = np.sqrt(np.square(y) + np.square(x))
        return np.exp(-r / H)

    def intgrl(x):
        if abs(x) > 1:
            integral = 2. * quad(integrand, 0, np.inf, args=(x, H))[0]
        else:
            lim = np.sqrt(1 - np.square(x))
            integral = quad(integrand, lim, np.inf, args=(x, H))[0]
        return integral

    xd = (np.arange(0, Ns*Rmax)+.5)/Ns

    Ich_r = np.vectorize(intgrl)(xd)
    Ich_tot = 2 * pi * sum(xd * Ich_r / Ns)  # int I 2 pi r dr

    if Ich_tot>0:
        Ich_r = Ich_r/Ich_tot

    if plotting:
        plt.plot(xd, Ich_r)
        plt.title('Chromosphere')
        plt.xlabel('r')
        plt.ylabel('F(r)')
        plt.show()

    return xd, Ich_r

def calc_disk_rprof(Rmax, Ns, plotting=1):
    nx = Ns*xmax*2
    xd = np.arange(0, nx)/Ns

    Fch_r = np.zeros(xd.shape)
    Fch_r[xd <= 1] = 1./(pi * Ns**2)

    if plotting:
        plt.plot(xd, Fch_r)
        plt.title('Chromosphere')
        plt.xlabel('r')
        plt.ylabel('F(r)')
        plt.show()

    return xd, Fch_r

# @njit(fastmath=True)
def calc_lc_from_rprof(xyp, Rp, Ir_fun, Ns, xpm, ypm, apm, vectorized=1):
    # xpm, ypm, apm: planet mask x,y and area
    # xyp: x,y position of planet in units of Rs
    # Ir_fun: function that returns brightness as function or r in units of Rs
    # Ns: stellar radius in pixels

    if vectorized:
        dpm = np.sqrt(np.square(xpm * Rp + np.atleast_2d(xyp[:, 0]).T) +
                      np.square(ypm * Rp + np.atleast_2d(xyp[:, 1]).T))
        fx_pl = np.sum(Ir_fun(dpm) * apm * (Ns * Rp) ** 2, axis=1)
    else:
        fx_pl = np.zeros(xyp.shape[0])
        for ipl in range(xyp.shape[0]):
            dpm = np.sqrt(np.square(xpm * Rp + xyp[ipl, 0]) + np.square(ypm * Rp + xyp[ipl, 1]))
            fx_pl[ipl] = np.sum(Ir_fun(dpm) * apm * (Ns * Rp) ** 2)

    fx_pl=fx_pl/Ns**2
    lc = 1 - fx_pl
    return lc

def tic():
    global tictocvar
    tictocvar=time.perf_counter()

def toc(str=''):
    global tictocvar
    print(f'{str}run time: {time.perf_counter()-tictocvar} sec')

def rms(a):
    return np.sqrt(np.mean(a**2))