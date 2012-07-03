#!/usr/bin/python
""" Prototype of "Show trajectories" in 2D in Python, no need to interact with C functions
1. opens parameter files and get the name of the files for _targets,
and numbers for rt_is.* and ptv_is.* files, as well as for images
2. loops through the sequence using the following set:
- opens the image
- opens the ptv_is.xxx and rt_is.xxx file (identical length and order)
- the 4 last columns in the rt_is. files are the number of lines in the respective 
_target files that give us the pixel coordinates to mark
- open the _targets file, last column is the number of the detected particle in the
rt_is. file - cross_validate and add a mark
- if in ptv_is. file this point is the beginning, i.e. the left is -1, then the color 
is red, if the right == -2, then the color is green, otherwise is a line
"""

import numpy as n
import read_parameters as read
import os
import matplotlib.pyplot as plt

n_cam = 4 # num of cameras


# get sequence path names, file names and values to build the loop
seq_names,first,last = read._read_sequence_par(n_cam) # 4 images, output is a list and two integers
# print seq_names, first, last

fig = plt.figure()


ax = []
for j in range(n_cam):
	axtmp = fig.add_subplot(2,2,j+1)
	axtmp.hold(True)
	ax.append(axtmp)

for i in range(first,last-1):
	rt_path = os.path.join('res','rt_is.'+str(i))
	ptv_path = os.path.join('res','ptv_is.'+str(i))
	rt_next_path = os.path.join('res','rt_is.'+str(i+1))
	ptv_next_path = os.path.join('res','ptv_is.'+str(i+1))
	
	# print rt_path, ptv_path
	rt = n.genfromtxt(rt_path,dtype=None,skip_header=1)
	ptv = n.genfromtxt(ptv_path,dtype=None,skip_header=1)
	rt_next = n.genfromtxt(rt_next_path,dtype=None,skip_header=1)
	ptv_next = n.genfromtxt(ptv_next_path,dtype=None,skip_header=1)
	
# 	if rt.size != ptv.size:
# 		print("rt_is and ptv_is have different number of particles")
#
	left = [row[0] for row in ptv]
	right = [row[1] for row in ptv]
	
	colors =[]
	for il in left:
		if il == -1:
			colors.append('r')
		else:
			colors.append('b')

	nn = [[row[4],row[5],row[6],row[7]] for row in rt]
	nn_next = [[row[4],row[5],row[6],row[7]] for row in rt_next]
	targets = []
	targets_next = []

	
	for j in range(n_cam):
		img_path = seq_names[j]+str(i)
		target_path = img_path+'_targets'
		# print img_path, target_path
		targets = n.genfromtxt(target_path,dtype=None,skip_header=1)
#
		img_path = seq_names[j]+str(i+1)
		target_path = img_path+'_targets'
		# print img_path, target_path
		targets_next = n.genfromtxt(target_path,dtype=None,skip_header=1)
		
		ax[j].imshow(plt.imread(img_path),origin='lower',alpha=0.5)
		for ir,row in enumerate(nn):
			if row[j] > -1: # there is a particle in this image
				x = targets[row[j]][1]
				y = targets[row[j]][2]
				if left[ir] == -1 and right[ir] != -2: # beginning, but not a single point
					ax[j].plot(x,y,'x',color=colors[ir]) # mark start by cross
					row_next = nn_next[right[ir]] 
					if row_next[j] > -1: # there is also a next point
						x1 = targets_next[row_next[j]][1]
						y1 = targets_next[row_next[j]][2]
						ax[j].plot([x,x1],[y,y1],'g-')
				elif left[ir] != -1 and right[ir] == -2: # end point, mark by blue cross
					ax[j].plot(x,y,'+',color=colors[ir]) # mark start by cross			
				elif right[ir] != -2: # all points that are trackable to the next step, get also lines
					row_next = nn_next[right[ir]] 
					if row_next[j] > -1: # there is also a next point
						x1 = targets_next[row_next[j]][1]
						y1 = targets_next[row_next[j]][2]
						ax[j].plot([x,x1],[y,y1],'g-')
			
		# ax[j].scatter(x,y,c=colors,s = 2)
		# ax[j].plot(x,y,'x-')
		
		
	plt.show()
			
			



