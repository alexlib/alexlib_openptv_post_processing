import pylab as p
import numpy as N
from tables import openFile

def plot_Lumley_triangle():
    ''' Plots Lumley triangle borders '''
    fig = p.figure()
    p.hold(True)
    x1 = N.ogrid[0:2/27.:50j]
    y1 = 3*(x1/2)**(2/3.)
    h1 = p.plot(x1,y1,'k')
    x2 = N.ogrid[-1/108.:0:50j]
    y2 = 3.*(x2/-2.)**(2/3.)
    h2 = p.plot(x2,y2,'k')
    # line([-1/108, 3*((1/108/2)^(2/3))],[2/27 3*((2/27/2)^(2/3))],'LineWidth',2,'Color','k');
    h3 = p.plot(N.ogrid[x2[0]:x1[-1:]:50j],N.ogrid[y2[0]:y1[-1:]:50j],'k')
    p.show()
    p.hold(False)
    return fig

def invariants(a):
    p = -N.trace(a) #sum(diag(a))
    q = 0.5 * (p**2 - N.trace(N.asmatrix(a)**2)) # sum(diag(power(a,2))))
    r = -N.linalg.det(a)
    return p,q,r
    
def anisotropy(x):
    # print 'x', x
    # print 'deviator', deviator(x)
    I,II,III = invariants(deviator(x))
    III = -III
    # print I, II, III
    return N.atleast_1d(I),N.atleast_1d(II), N.atleast_1d(III)
    
def deviator(x):
    x = x.astype('f')
    return x/N.trace(x) - N.eye(3)/3.0 # x/sum(diag(x)) - eye(3)/3.0
    
    
fig = plot_Lumley_triangle()
p.hold(True)
a = N.array([[1,0,0],[0,0,0],[0,0,0]]).astype('f')
I,II,III = anisotropy(a);
p.plot(III,-II,'o',label='1')
a = N.array([[1,0,0],[0,1,0],[0,0,0]]).astype('f')
I,II,III = anisotropy(a);
p.plot(III,-II,'o',label='2')
a = N.array([[1,0,0],[0,1,0],[0,0,1]]).astype('f')
I,II,III = anisotropy(a);
p.plot(III,-II,'o',label='3')
a = N.array([[1,1,0],[1,0,0],[0,0,0]]).astype('f')
I,II,III = anisotropy(a);
p.plot(III,-II,'o',label='4')
a = N.array([[1,1,0],[1,1,0],[0,0,0]]).astype('f')
I,II,III = anisotropy(a);
p.plot(III,-II,'o',label='5')
a = N.array([[1,1,0],[1,1,0],[0,0,1]]).astype('f')
I,II,III = anisotropy(a);
p.plot(III,-II,'o',label='6')
a = N.array([[1,1,1],[1,0,0],[1,0,0]]).astype('f')
I,II,III = anisotropy(a);
p.plot(III,-II,'o',label='7')
a = N.array([[1,1,1],[1,1,0],[1,0,0]]).astype('f')
I,II,III = anisotropy(a);
p.plot(III,-II,'o',label='8')

a = N.array([[1,1,1],[1,1,0],[1,0,1]]).astype('f')
I,II,III = anisotropy(a);
p.plot(III,-II,'.',label='9')
a = N.array([[1,1,1],[1,0,1],[1,1,0]]).astype('f')
I,II,III = anisotropy(a);
p.plot(III,-II,'.',label='10')
a = N.array([[1,1,1],[1,1,1],[1,1,0]]).astype('f')
I,II,III = anisotropy(a);
p.plot(III,-II,'.',label='11')

a = N.array([[1,1,1],[1,1,1],[1,1,1]]).astype('f')
I,II,III = anisotropy(a);
p.plot(III,-II,'.',label='12')
