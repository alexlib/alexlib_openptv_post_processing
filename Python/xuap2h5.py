# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Documents and Settings\user\.spyder\.temp.py
"""
import sys
import getopt
import glob
import os

def main(argv=None):
    if argv is None:
        argv = sys.argv
        print argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.error, msg:
             raise Usage(msg)
        # more code, unchanged
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2
        
    readxuapfiles(argv)
    
def readxuapfiles(directory):
    print "hello"
    lst=glob.glob(os.path.join(directory,"xuap.*"))
    print lst
    data = [[] for i in range(len(lst))]
    for i,name in enumerate(lst):
	tmp = n.loadtxt(name)
	data[i] = n.rec.fromrecords(tmp.tolist(),names='l,r,xr,yr,zr,xf,yf,vf,u,v,w,ax,ay,az,s2n')
    return data    
        
        
if __name__ == "__main__":
    sys.exit(main())