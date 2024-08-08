## TACHELES

  TACHELES (Transits Across CHromosphEricaLly/coronally activE Stars) is a model for Exoplanet transits which includes, beside a photosphere, a numerical implementation of the chromosphere/corona, i.e. the optically thin outer layer of a star.

#Installation
------------
  To install the required packes, run

	'''
   	pip3 install PyAstronomy astropy math numba numpy photutils scipy batman csv matplotlib os time
   	'''

  To install the latest release via pip from PyPI use
  
    '''
    pip install TACHELES
    '''
    
  to install the current version from github. 
  
  Alternatively, download the source and use
  
    '''
    python -m pip install .
    '''
    
#Parameters
-----------
Aside from the standard parameters of an Exoplanet system (Rp, period, t0, a, e, inc, w, Omega, u1, u2) implemented via BATMAN, TACHELES requires several additional parameters for the modeling of the outer stellar atmosphere. 

Since the integrals necessary to calculate the line-of-sight flux are analytically unsolvable, TACHELES utilizes a numerical estimation. 
In order to do this, TACHELES first performs a grid calculation, and then interpolates on that grid. It requires input on the grid parameters:

dlogH: fineness of the grid in logarithmic scale height
logHup: upper limit of the grid in logH
logHlow: lower limit of the grid in logH
Rmax: extend of the image in stellar radii
npl: number of planet pixels
osf: oversample factor
npl: stellar radius in pixels

For all of these values, the defaults are the ones used in Perdelwitz et al. (2024).

#Quickstart
-----------
In order to compute a TACHELES light curve using the system parameters of HAT-P-18 b, 

'''
from TACHELES.TACHELES import *
import numpy as np

for i in range(-3,1):
	e = 0.
	w = 90.
	Omega = 90.
	a = 16.52
	inc = 88.9
	t0 = 0.
	period = 5.5080232
	logH=-1.5
	logBr=-0.5
	u1=0.7
	u2=0.2
	Rp=0.139
	t=np.linspace(-0.1,0.1,10000)

	lc=TACHELES_lc(t, Rp, period, t0, a, e, inc, w, Omega, logH, float(i), u1, u2)

	plt.plot(t,lc,label='logBr='+str(i))
plt.legend(loc='lower left')
plt.xlabel('time [d]')
plt.ylabel('flux [norm]')
plt.show()
'''

TACHELES first computes a chromospheric grid (which may take a while depending on the grid parameters) and uses it to produce a light curve for each of four different brightness ratios, which are then plotted.
The grid for any given set of parameters is calculated only once and saved in the same directory. TACHELES automatically checks whether a suitable grid is already available in the directory, and uses it instead of recalculating the grid.


#Acknowledgements and Citation
-----------------------------
If you use TACHELES for your publication, please cite the following paper:
Perdelwitz et al. 2024 ...

Licensing
---------

  Where not stated otherwise, TACHELES is released under the
  MIT license.
