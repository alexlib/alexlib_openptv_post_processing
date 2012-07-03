from numpy import *
from matplotlib.mlab import find
from pylab import *

B = loadtxt("77000000.txt")
x = B[:,0]
y = B[:,1]
u = B[:,2]
v = B[:,3]

x = x[x!=0]
ux = unique(x)
minx = min(ux)
maxx = max(ux)

y = y[y!=0]
uy = unique(y)
miny= min(uy)
maxy = max(uy)



xx,yy = meshgrid(linspace(minx,maxx,ux.shape[0]),linspace(miny,maxy,uy.shape[0]))
rows,cols = xx.shape

uu = zeros_like(xx)
vv = zeros_like(yy)

for i,n in enumerate(x):
    m,n = unravel_index(find(logical_and((xx==x[i]),(yy==y[i]))),(rows,cols))
    uu[m,n] = u[i]
    vv[m,n] = v[i]



ux,uy = gradient(uu)
vx,vy = gradient(vv)
vort = vx - uy


M = zeros(uu.shape, dtype='bool')
M[xx<900] = True
M[xx>1550] = True
M[abs(vort)>2]=True
M[yy<200] = True

uu = ma.masked_array(uu,mask=M)
vv = ma.masked_array(vv,mask=M)
vort = ma.masked_array(vort,mask=M)


figure()
im = imread('000077RA.TIF')
imshow(im,origin='lower',cmap=cm.gray)
hold('on')
Q = quiver(xx,yy,uu,vv,vort,pivot='middle',scale=400,width=.0015)
clim(-1.5,1.5) 
colorbar()
# clim(-.1,.1)
  