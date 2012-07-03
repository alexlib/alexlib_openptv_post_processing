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
    
def invariants(x):
    p = -N.trace(x) #sum(diag(a))
    q = 0.5 * (N.power(p,2) - N.trace(N.asmatrix(x)**2)) # sum(diag(power(a,2))))
    r = -N.linalg.det(x)
    # print p,q,r
    return p,q,r


def deviator(x):
    # x = x.astype('f')
    return x/N.trace(x) - N.eye(3)/3.0 # x/sum(diag(x)) - eye(3)/3.0   
    
def anisotropy(x):
    # print 'x', x
    # print 'deviator', deviator(x)
    I,II,III = invariants(deviator(x))
    III = -III
    return N.atleast_1d(I),N.atleast_1d(II), N.atleast_1d(III)
    



   

# def calculate_save_anisotropy_invariants()

# for k in range(dataNumVec):

    # j = dataNumVec(k); 

    # load(['meanCorField_',int2str(j)]);
    # [II,III,Au] = deal(zeros(length(meaCor.uu),1));
    # for i=1:length(meaCor.uu),
        # [II(i),III(i)] = anisotropy([meaCor.uu(i),meaCor.uv(i),meaCor.uw(i);
            # meaCor.uv(i),meaCor.vv(i),meaCor.vw(i);
            # meaCor.uw(i),meaCor.vw(i),meaCor.ww(i)]);
        # % scatter(iii,-ii,'r.');
        # Au(i) = (III(i)/2)/((-II(i)/3)^(3/2));
        # F(i) = 9*II(i) + 27*III(i)+1; 
    # end
    # save(['anisotropyCorField_',int2str(j)],'II','III','Au','F');


    




h5file = openFile('run4_trajPoint.h5',"r")
data = h5file.root.trajPoint

fig = plot_Lumley_triangle()
p.hold(True)
for row in data[1000:1200]:    
    u = row['u']
    v = row['v']
    w = row['w']
    # print u, v, w
    #Rs = N.array([[u*u, u*v, u*w],[v*u, v*v, v*w],[w*u, w*v, w*w]]).astype('f')
    Rs = N.array([[u*u, u*v, u*w],[v*u, v*v, v*w],[w*u, w*v, w*w]])
    # print Rs/N.trace(Rs)
    # print deviator(Rs)
    I,II,III = anisotropy(Rs)
    # print I, II, III
    p.plot(III,-II,'o')
    # del I, II, III
    
    


def test_lumley_triangle():
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

h5file.close()

