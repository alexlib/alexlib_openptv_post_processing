"""
RQ plot
"""
from tables import openFile
from numpy import r_, c_, mean, trace, matrix, histogramdd, log, sum, zeros
from pylab import figure, plot, contourf, colorbar, hold, xlabel, ylabel, axis, show 

h5file = openFile('../080825/run6_080825.h5',"r")
data = h5file.root.trajPoint

# Empirics
mean_s = 2.8

Q = zeros((len(data),1),'f')
R = zeros((len(data),1),'f')

for k,row in enumerate(data[::100]):
    ux = row['s11']
    uy = row['s12']-0.5*row['w3']
    uz = row['s13']+0.5*row['w2']
    vx = row['s12']+0.5*row['w3'] 
    vy = row['s22']
    vz = row['s23']-0.5*row['w1']
    wx = row['s13']-0.5*row['w2']
    wy = row['s23']+0.5*row['w1']
    wz = row['s33']    
    
    uM = r_[ux,uy,uz, vx, vy,vz,wx,wy,wz].reshape(3,3)
    tmpQ = -(1/2.)*(trace(matrix(uM)**2))/mean_s**2
    tmpR = -(1/3.)*(trace(matrix(uM)**3))/mean_s**3
    if abs(tmpQ) <= 100 and abs(tmpR) <= 100 and sum(uM**2) <= 100:
        Q[k] = tmpQ
        R[k] = tmpR
    
    

figure()
hold('True')
jRQ,xedges = histogramdd(c_[Q,R],200)
contourf(log(jRQ+1)/log(10.),50)

xlabel('R')
ylabel('Q');
colorbar()
axis([70, 150, 40, 150])
show()    
