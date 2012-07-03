#!/usr/bin/env python
"""
Author: Alex Liberzon
Date:   22-Mar-2011

Renames TIFF files acquired with Streams5 suitable for Insight.
Usage: 
	python rename.py ./scene21
"""

import os
import sys
# from optparse import OptionParser


import os
directory = sys.argv[1] # './scene21'
# import pdb; pdb.set_trace()
if not os.path.isdir(directory):
    print "error in directory"

filelist = os.listdir(directory)
filelist.sort()
# using only pairs, sub sampling at half-frequency
filelistA = filelist[0:-1:2]
filelistB = filelist[1::2]

##  if both lists are duplicate
# filelistA = filelist[:-1]
# filelistB = filelist[1:]

replacable = '_1-DVR Express CLFC_'

for filename in filelistA:
    new = filename.replace(replacable,'_')
    new = new.upper().replace('.TIF','.LA.TIF')
    old = os.path.join(os.path.abspath(directory),filename)
    os.rename(old,os.path.join(os.path.abspath(directory),new))

for filename in filelistB:
    new = filename.replace(replacable,'_')
    new = new.upper().replace('.TIF','.LB.TIF')
    i = new.rfind('_')+1
    j = new[i:].find('.')
    new = new.replace(new[i:i+j], '%04d' % (int(new[i:i+j]) - 1))
    old = os.path.join(os.path.abspath(directory),filename)
    os.rename(old,os.path.join(os.path.abspath(directory),new))


