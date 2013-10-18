# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
from glob import glob
import numpy as np
from matplotlib.mlab import csv2rec
from mpl_toolkits.mplot3d.axes3d import Axes3D

# <codecell>

directory = '/Users/alex/Documents/PTV/test/res/'

# <codecell>

list_ptv_is_files = glob(os.path.join(directory,'ptv_is.*'))
# print list_ptv_is_files

# <codecell>

# reading the ptv_is.* files
frames = []
for counter, ptv_is_file in enumerate(list_ptv_is_files):
    frame = csv2rec(ptv_is_file,skiprows=1,delimiter=' ',names=['p','n','x','y','z'])
    frame = rec_append_fields(frame,['t','id'],[np.zeros_like(frame.x)+counter,np.zeros_like(frame.x)-999],dtypes=[np.int,np.int])
    frames.append(frame)
    

# <codecell>

# adding trajectory id = linking
id = 0
for i,f in enumerate(frames):
    for j,l in enumerate(f):
        if l.p == -1 and l.n != -2:
            l.id = id
            id += 1
        elif l.p != -1:
            l.id = frames[i-1].id[l.p]
        

# <codecell>

for i,f in enumerate(frames):
    ind = f.id == -999
    frames[i] = f[~ind]

# <codecell>

last_traj = max(frames[-1].id)
traj = [[] for k in range(last_traj+1)]
for f in frames:
    for p in f:
        traj[p.id].append(p)

# <codecell>

fig = figure(figsize=(12,8))
ax = fig.add_subplot(1,1,1, projection='3d')
for t in traj:
    x = [p.x for p in t]
    y = [p.y for p in t]
    z = [p.z for p in t]
    ax.plot(x,y,z)

