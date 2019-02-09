import sys
import numpy as np
from mfd import sca
from matplotlib import pyplot as pl

n = int(sys.argv[1])
xr = np.linspace(0, 2, n)
yr = xr[:int(n*.66)]
w = abs(xr[1]-xr[0])
x, y = np.meshgrid(xr, yr)
r = np.sqrt(x*x+y*y)
z = np.exp(-r*r)
z[r > 1.75] = np.nan
a = sca(z, w)
e = 100 * (2 * a - r) / r
e = np.ma.masked_invalid(e)

pl.figure(1, (13.66, 7.68))
pl.title('Gaussian hill')
pl.pcolormesh(x, y, e,
        vmin = -100, vmax = 100,
        cmap = pl.cm.seismic)
cb = pl.colorbar()
cb.set_label(r'Rel. deviations [%]')
pl.axes().set_aspect('equal')
pl.savefig('%s.png' % sys.argv[0][:-3])
