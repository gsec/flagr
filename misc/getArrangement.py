import numpy as np
import pylab
import matplotlib.pyplot as plt

N = 10     # number-"radius" of the arrangement

def getArrangement(geo='quad',bounding='rect', offset=(0,0), plot=0):
    #builds matrix of particles location 
    #1st argument: quad or hex, 2nd: rect or circle
    if geo=='hex':
        Nmax = int(N*np.sqrt(3))
        #assert Nmax % 2 == 0, "Nmax is odd!" 
        nums = np.arange((-Nmax/2),(Nmax/2+1),1.0)
        mx, ny = np.meshgrid(nums,nums)
        mx = mx + 0.5*ny + offset[0]
        ny = ny * np.sqrt(3)/2. + offset[1]
    else:
        nums = np.arange((-N/2),(N/2+1),1.0)
        mx, ny = np.meshgrid(nums,nums)
    if bounding=='circle':
        ind = np.nonzero(mx**2 + ny**2 <= (N/2)**2)
    else:
        ind = np.nonzero(mx**2 + ny**2 <= np.inf)
    m = mx[ind]
    n = ny[ind]
    #m = np.delete(m,(len(m)-1)/2)
    #n = np.delete(n,(len(n)-1)/2)
    if plot==1:
        plt.axis("equal")
        plt.plot(m,n,'o')
        plt.show()
    return m, n
    
a = np.array((0,0))
b = (1./2)*np.array((1., 1./np.sqrt(3)))
c = 2*b

m1,n1 = getArrangement('hex','circle', offset=a)
m2,n2 = getArrangement('hex','circle', offset=b)
m3,n3 = getArrangement('hex','circle', offset=c)

plt.axis("equal")
plt.plot(m1,n1,'o')
plt.plot(m2,n2,'o')
plt.plot(m3,n3,'o')
plt.show()
