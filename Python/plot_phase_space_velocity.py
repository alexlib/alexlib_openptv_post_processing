import pylab as p
import numpy as N
from tables import openFile

h5file = openFile('run4_trajPoint.h5',"r")
data = h5file.root.trajPoint

# fig = plot_Lumley_triangle()

p.figure()
p.hold(True)

# for row in data[1000:1200]:    
u = data[:]['u']
v = data[:]['v']
w = data[:]['w']
p.subplot(311)
p.plot(u,v,'.')
p.subplot(312)
p.plot(u,w,'.')
p.subplot(313)
p.plot(v,w,'.')
    
    