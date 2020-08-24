from numpy import *
from pylab import *

figure(1)
clf()

x = loadtxt('out.dat')
plot(x[:,0],x[:,1],'k-')

y = loadtxt('bunny10')
plot(y[:,0],y[:,1],'rx')
show()
