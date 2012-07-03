"""
build_trajectories_from_xuap(data) is responsible for the book keeping of the trajectories, 
using the prev,next fields in the frames,creating new set of data, in the form of trajectories
see also ptv_is_to_traj for _build_trajectories _from ptv_is
"""
import numpy as np

def build_trajectories_from_xuap(data):
    # first frame, initialization
    trajid = 0 # running trajid counter
    data[0].trajid = np.nan*np.empty_like(data[0].xr)
    data[0].t = data[0].t*np.ones_like(data[0].xr)
    ind = data[0].next > -1
    if data[0].trajid.size > 1:
    	data[0].trajid[ind.nonzero()] = range(trajid,trajid+ind.size)
    else:
    	data[0].trajid = trajid
    	
    trajid = trajid+ind.size
    for i,frame in enumerate(data[1:]):
        frame.trajid = np.nan*np.empty_like(frame.xr)
        old = frame.prev > 0
        frame.trajid[old] = data[i].trajid[frame.prev[old]] # notice i runs from zero
        ind = np.logical_and(frame.prev < 1, frame.next > -1)
        if data[0].trajid.size > 1:
            frame.trajid[ind.nonzero()] = range(trajid,trajid+ind.size)
        else:
            frame.traji = trajid
            
        trajid = trajid+ind.size
        frame.t = frame.t*np.ones_like(frame.xr)
        
    return(data)
    
if __name__ == "__main__":
	npz = np.load('markxuapdata.npz')
	data = npz['arr_0']
	newdata = build_trajectories_from_xuap(data)
