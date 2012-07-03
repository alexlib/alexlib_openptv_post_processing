import numpy as np
import pylab as pl

N = 1000
n = 10
np.random.seed(3)#use always the same seed
x, y = np.random.randn(2, N)/10 +0.5
X, Y = np.mgrid[0:1:n*1j, 0:1:n*1j]

xfloor = X[:,0][np.floor(n*x).astype(int)]
yfloor = Y[0][np.floor(n*y).astype(int)]
z = xfloor + n*yfloor
Z = X + n*Y
histo = np.histogram(z.ravel(), bins=np.r_[Z.T.ravel(),2*n**2])

pl.pcolor(X-1./(2*n), Y-1./(2*n), histo[0].reshape((n,n))) #shifted to
# have centered bins
pl.scatter(x, y)
pl.show()