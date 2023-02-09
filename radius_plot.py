import numpy as np
import matplotlib.pyplot as plt
from sys import exit


Ndata = 500
Ravg = []
for i in range(Ndata):
    data = np.loadtxt("data/drops_{}.dat".format(i))[:,2]
    R = np.asarray(data[ data > 0] )
    Ravg.append(np.mean(R))
   
plt.plot(Ravg) 

xmin = 1.0
xmax = 2*Ndata
x = np.linspace(xmin,xmax, 1000)
y = (0.5) * x**(1.0/3)
plt.plot(x,y)

plt.xlim([xmin, xmax])
#plt.ylim([0.5, 1])
plt.xscale("log")
plt.yscale("log")
plt.show()
