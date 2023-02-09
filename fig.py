import numpy as np
import matplotlib.pyplot as plt
from sys import exit
Ndata = 500
plt.subplot(2,1,1)
t_list = [0, 250, 500]
for t in t_list:
    data = np.loadtxt("data/drops_{}.dat".format(t))[:,2]
    R = np.asarray(data[ data > 0] )

    h, bedges = np.histogram(R, density=True, bins = 10) 

    db = bedges[1] - bedges[0]
    bedges += db /2
    bins = bedges[:-1]
    plt.plot(bins, h, label=t)

plt.legend()
plt.subplot(2,1,2)



Ravg = []
for t in range(Ndata):

    data = np.loadtxt("data/drops_{}.dat".format(t))[:,2]
    R = np.asarray(data[ data > 0] )

    Ravg.append(np.mean(R))
   
plt.plot(Ravg)

x = np.linspace(1, Ndata,Ndata)
y = 1 * x**(1.0/3.0)
plt.plot(x,y)

plt.yscale("log")
plt.xscale("log")

plt.show()
