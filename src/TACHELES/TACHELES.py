from PyAstronomy import pyasl
from astropy.io import fits
from TACHELES.chrom_lib import *
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
import numpy as np



def TACHELES_lc(t, Rp, period, t0, a, e, inc, w, Omega, logH, logBr, u1, u2, off=1, dlogH=0.05, logHup=0.05, logHlow=-3.05, Rmax=8, Ns=1000, npl=100, osf=1):
	# t: time grid (np.array)
	# Rp: planetary radius in units of stellar radii
	# t0: transit midpoint
	# period, a, e, inc, w, Omega: Keplerian elements, as implemented in BATMAN
	# logH: logarithmic scale height of the chromosphere/corona, in units of stellar radii
	# logBr: brightness ratio of the chromospheric/coronal disk with regard to the photospheric disk
	# u1, u2: quadratic limb darkening coefficients for the photosphere
	# off: multiplicative offset of the light curve
	# dlogH: fineness of the logarithmic scale height of the grid
	# logHup: upper limit for the logarithmic scale height of the grid
	# logHlow: lower limit for the logarithmic scale height of the grid
	# Rmax: Extend of the 2D image and light curve in units of stellar radii
	# Ns: stellar radius in pixels
	# npl: number of pixels in the planetary disk
	# osf: oversampling factor

	params = batman.TransitParams()
	params.t0, params.per, params.rp, params.a, params.inc, params.ecc, params.w, params.u, params.limb_dark = t0, period, Rp, a, inc, e, w, [u1, u2], "quadratic"
	m = batman.TransitModel(params, t)
	flux_curve = m.light_curve(params)
	lc=(np.power(10,logBr)*lc_chrom_calc(t, Rp, logH, dlogH, logHup, logHlow, Rmax, Ns, period, a, inc, e, t0, w, Omega, npl, osf)+flux_curve)/(1.+np.power(10,logBr))
	total_flux = np.nan_to_num(lc, nan=1)
	lc=total_flux*off
	return lc

