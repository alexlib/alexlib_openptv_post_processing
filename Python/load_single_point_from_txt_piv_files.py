from scipy import * 
from pylab import *
import os
import glob
"""
Test the fft routine. Add signals, and multiply signals. 
"""

# load data:

# a = loadtxt('scene23_000000.txt')
# find((a[:,0] == 1008) & (a[:,1] == 784)) # or 736
# the answer is 1154, or 1388

# 
u = []; v = []
directory = '/Users/alex/Desktop/scene23'
for i in glob.glob(os.path.join(directory,'*.txt')):
    _ = loadtxt(i,skiprows=1388)
    u.append(_[0,2])
    v.append(_[0,3])
    
# savetxt('u1154.txt',u);
# savetxt('v1154.txt',v);

u = asarray(u)
v = asarray(v)

def nextpow2(n):
    return 2**(ceil(log2(n))-1)


npts = nextpow2(len(u)) # 512            #Use some power of 2
t=linspace(0,1,npts+1)     # Use 2^N + 1 
dt = (t[-1]-t[0])/(len(t)-1)    # Maximum frequency is 1/2dt ?
fmax = 1/(2*dt) 
# f1 = 80 
# f2 = 90 
#sig = 1 + sin(2*pi*f1*t) + 1 + sin(2*pi*f2*t)  # sum of signals 
# sig=(1+sin(2*pi*f1*t))*(1+sin(2*pi*f2*t))      # product of signals
sig = u[:npts+1] - mean(u)

figure(1)
plot(t,sig);xlabel('Time');title('Signal'); show()

ft = fft(sig,n=npts) 
mgft=abs(ft)             #Get magnitude of fft
df = fmax/float(npts/2)
f=linspace(0,fmax,npts/2+1)
print 'fmax = ',fmax,' df = ',df,' ','\n 1st freqs = ',f[0:5]

figure(2)
plot(f,mgft[0:npts/2+1]);title('Fast Fourier Transform Magnitude')
xlabel('frequency')
ylabel('fft magnitude')
show()
