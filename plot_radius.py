import numpy as np
import matplotlib.pyplot as plt
from sys import exit



data = np.loadtxt("mean_radius.dat");
t = data[:,0]
radius = data[:,1]
number = data[:,2]


y = 0.1 * t**(1.0/3.0)

plt.subplot(2,1,1)
plt.plot(t, radius)
plt.plot(t,y)

plt.xscale("log")
plt.yscale("log")
plt.ylabel(r"$<R>$", rotation=0)

plt.subplot(2,1,2)
plt.plot(t, number)
plt.yscale("log")
plt.xscale("log")

plt.xlabel("time")
plt.ylabel(r"$N$", rotation=0)

plt.tight_layout()
plt.savefig("averageR.pdf")
#plt.show()
