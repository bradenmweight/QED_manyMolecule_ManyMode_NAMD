import numpy as np
from matplotlib import pyplot as plt

NDirs = 100
NSteps = 1000

TEMP = np.zeros(( NDirs, NSteps ))
for d in range( NDirs ):
    TEMP[d,:] = np.loadtxt(f"traj-{d}/MD_OUTPUT/Temperature.dat")[:NSteps,1]

print("Average Temperature (Expected) = %2.4f K (%2.0f K)" % (np.average(TEMP), 300) )
print("STD     Temperature = %2.4f K" % (np.std(TEMP)) )

for d in range(NDirs):
    plt.plot( np.arange(1,NSteps+1), TEMP[d,:], alpha=0.25, c="black", lw=2 )
plt.plot( np.arange(1,NSteps+1), np.average(TEMP[:,:],axis=0), c="red", lw=4, label="Average" )

plt.legend()
plt.xlabel("MD Step", fontsize=18)
plt.ylabel("Temperature (K)", fontsize=18)
plt.savefig("temperature.jpg",dpi=600)

