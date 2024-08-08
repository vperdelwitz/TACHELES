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