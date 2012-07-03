from tables import openFile
import numpy as np
from pylab import figure, plot, show, hold

from scipy import fftpack
import math

def ExtractFFt(time,data):

    N = len(time)
    P1 = (time[N-1]-time[0])
    n1 = 1.0/P1;

    dft = fftpack.fft(data)
    # print dft
    x = []
    y = []

    for i in range(0,N/2):
       y.append(math.sqrt(np.real(dft[i])**2 + np.imag(dft[i])**2)*2.0/N )
       x.append(n1*i)

    return(x,y)
   

h5file = openFile('run2_080827.h5',"r")
# h5file = openFile('run3_080825.h5',"r")


data = h5file.root.trajPoint

u = data[::10]['u']
v = data[::10]['v']
w = data[::10]['w']
t = data[::10]['t']

h5file.close()

uniqueT = np.unique(t)
umean = vmean = wmean = np.empty_like(uniqueT).astype('f')

for k,i in enumerate(uniqueT):
    umean[k] = np.mean(u[t==i])
    vmean[k] = np.mean(v[t==i])
    wmean[k] = np.mean(w[t==i])

# figure()
# plot(t,u,'b.')
# hold('True')

# show()

# figure()
# plot(t,v,'b.')
# hold('True')
# plot(uniqueT,vmean,'r.')
# show()

# figure()
# plot(t,w,'b.')
# hold('True')
# plot(uniqueT,wmean,'k.')
# show()

figure()
hold('True')
plot(uniqueT,umean,'g.')
plot(uniqueT,vmean,'r.')
plot(uniqueT,wmean,'k.')
show()


f,sp = ExtractFFt(uniqueT,umean)

figure()
plot(f[:],sp[:])
hold('True')
f,sp = ExtractFFt(uniqueT,vmean)
plot(f[:],sp[:])
f,sp = ExtractFFt(uniqueT,wmean)
plot(f[:],sp[:])
show()
 


