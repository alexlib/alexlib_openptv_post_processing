#!/usr/bin/python
""" analyse the crowd tracking data
"""
import os
import numpy as np
import matplotlib.pyplot as plt
os.chdir('/Users/alex/Documents/PTV/ptv_postproc/Python/')

from crowd_ptv_is_to_traj import *
data = np.load('/Users/alex/Desktop/crowd_tracking/ptv_is.npz')

traj = data['traj']
data = data['data']

density = []
mean_abs_u = []
mean_abs_v = []
mean_var_u = []
mean_var_v = []

for f in data:
    density.append(len(f))
    
    # local mean velocity
    u = 0.0
    v = 0.0    
    for p in f:
        u += np.abs(p.u) 
        v += np.abs(p.v)
        
    # mean absolute amplitude    
    u = u/len(f)
    v = v/len(f)
    
    # now the local fluctuations (some kind of a stress)
    uf = 0.0
    vf = 0.0
    for p in f:
        uf += (np.abs(p.u) - u)**2 #sort of local variance
        vf += (np.abs(p.v) - v)**2
    
    mean_var_u.append(uf/len(f))
    mean_var_v.append(vf/len(f))
         
    mean_abs_u.append(u/len(f))
    mean_abs_v.append(v/len(f))

plt.figure()
plt.plot(density,mean_abs_u,'o--')
plt.plot(density,mean_abs_v,'s:')
plt.title('Local mean amplitude of U,V vs density')
plt.xlabel('Density [#]')
plt.ylabel('Velocity amplitude [m/s]')
plt.legend(('U','V'))
plt.axis([0,6, -.1, 2])
plt.show()

plt.figure()
plt.plot(density,mean_var_u,'o')
plt.plot(density,mean_var_v,'s')
plt.title('Local mean variance vs density')
plt.xlabel('Density [#]')
plt.ylabel('Velocity amplitude variance [m/s]')
plt.legend(('U','V'))
plt.show()



# how to get mean velocity field? 
# we need for each x region a value of velocity
# some sort of Eulerian interpolation

x = []
y = []
u = []
v = []

for f in data:
    for p in f:
        x.append(p.x)
        y.append(p.y)
        u.append(p.u)
        v.append(p.v)

minx = np.min(x)
maxx = np.max(x)
miny = np.min(y)
maxy = np.max(y)

# let's try griddata
from scipy.interpolate import griddata
[X,Y] = np.meshgrid(np.linspace(minx,maxx,100),np.linspace(miny,maxy,100))
U = griddata(x,y,u,X,Y)
plt.figure()
plt.pcolor(X, Y, U, cmap=cm.jet)

# let's try ensemble average

U = np.zeros((100,100))
for f in data: # for each frame
    xt = []
    yt = []
    ut = []
    vt = []
    for p in f:
        xt.append(p.x)
        yt.append(p.y)
        ut.append(p.u)
        vt.append(p.v)
    
    xt = np.asarray(xt)
    yt = np.asarray(yt)
    ut = np.asarray(ut)
    vt = np.asarray(vt)
        
    # now fit to the grid
    tmp = griddata((xt,yt),ut,(X,Y),method='nearest')
    U += tmp

U = U/len(data) # from sum to average
U = U - mean(U.flatten()) # centralize around zero

fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolor(X, Y, U, cmap=cm.jet)
plt.axis([-1.04,1.34,-.5,.4])
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.colorbar(orientation='horizontal')
plt.title('Average velocity map')
ax.set_aspect('equal')
fig.show()

    
    
    


