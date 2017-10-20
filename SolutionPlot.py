###Plot solution of ICG model ###

import scipy as SP
import matplotlib as MP
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import numpy as np
import math as m

#Load solution (from the model, computed with Model,c)
Sol = SP.loadtxt('./Res/data.txt')
Sol = SP.array(Sol)

time = Sol[:,0]/60.0
Nblood = Sol[:,1]
Ns = Sol[:,2]
Nh = Sol[:,3]
Nbc = Sol[:,4]
Liver = Ns + Nh + Nbc

Fig1 = plt.figure()
ax_1 = Fig1.add_subplot(311)
ax_2 = Fig1.add_subplot(312)
ax_3 = Fig1.add_subplot(313)

ax_1.plot(time, Nblood, '-r',linewidth = 1.5, label = 'Blood circulation')

ax_2.plot(time, Ns, '-m', linewidth = 1.5, label = 'Sinusoids')
ax_2.plot(time, Nh, '-b', linewidth = 1.5, label = 'Hepatocytes')
ax_2.plot(time, Nbc, '-g', linewidth = 1.5, label = 'Bile canaliculi')

ax_3.plot(time, Liver, '-g', linewidth = 1.5, label = 'Liver')

ax_1.set_xlabel("time in min");
ax_2.set_xlabel("time in min");
ax_3.set_xlabel("time in min");

ax_1.set_ylabel("ICG qty (AU.ml)");
ax_2.set_ylabel("ICG qty (AU.ml)");
ax_3.set_ylabel("ICG qty (AU.ml)");

ax_1.legend(loc = 'best')
ax_2.legend(loc = 'best')
ax_3.legend(loc = 'best')


# Comparison to observation
#Load Observation (from El-Desoky et al. paper)
Obsdata = SP.loadtxt('./Obs/ElDes_Control_qty.txt')
obsTime = SP.array(Obsdata[:,0])/60.0
obsLiv = SP.array(2**Obsdata[:,1])


Fig2 = plt.figure()
ax_4 = Fig2.add_subplot(111)

ax_4.plot(obsTime,obsLiv, '.', color = '0.75', label = "Observation")
ax_4.plot(time,Liver, '-r', linewidth=1.5, label = "Model")

ax_4.legend(loc = 'best')
ax_4.set_xlabel("time in s");
ax_4.set_ylabel("ICG in liver qty (AU.ml)");

plt.show()
