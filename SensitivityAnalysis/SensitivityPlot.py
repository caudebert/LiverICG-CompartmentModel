"""
============================================================================
This file is part of the code ICG4pharma
ICG 4-compartment liver pharmacokinetics model
Copyright (C) 2017, INRIA
*******************************************************
******* Author : Chloe Audebert and Irene E. Vignon Clementel
******* Last mod : 5/01/2016
*******************************************************
******* Description ***********
Model is considered to fit El-Desoky et al. 1999  data of ICG fluorescence in liver rabbits
4 compartments are considered, the liver with 3 compartments : Sinusoid, hepatocytes and bile canaliculi 
and the rest of the blood circulation
Concentrations are in arbitrary units
Initial value in the blood is assumed non-zero (just after injection)

The sensitivity equation are solved in Model_Sensitivity.c
Observation : Liver amount
Parameter of interest : Ksh Qhb S and Fb

This python script plot the traditional and generalized sensitivity function (from the output file of Model_Sensitivity.c)

ICG4Pharma is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.
ICG4Pharma is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
more details.
You should have received a copy of the GNU Lesser General Public License
along with ICG4Pharma. If not, see http://www.gnu.org/licenses/.
 =============================================================================
 """
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from matplotlib import rc

raw_sens_liv = sp.loadtxt("./Sens_L.txt")
raw_sens_liv = sp.array(raw_sens_liv)
scale_liv = sp.loadtxt("./D_L.txt")
sl1 = raw_sens_liv[:,0]*scale_liv[0]; #Ksh
sl2 = raw_sens_liv[:,1]*scale_liv[1]; #Qhb
sl3 = raw_sens_liv[:,2]*scale_liv[2]; #Fbc
sl4 = raw_sens_liv[:,3]*scale_liv[3]; #S

#solution
Sol = sp.loadtxt("./data.txt");
time = Sol[:,0]/60.0
liv = Sol[:,5]

#plot traditional sensitivity
Fig = plt.figure(figsize = [12.,6.])
a = Fig.add_subplot(121)
a2 = Fig.add_subplot(122)
a.plot(time, liv, linewidth = 1.5, label = 'liver ICG  amount')
a.legend
a.set_ylabel('ICG amount (AU.ml) ')
a.set_xlabel('time (s)')

a2.plot(time, sl1, 'k', linewidth = 1.5, label = 'Ksh')
a2.plot(time, sl2, 'r', linewidth = 1.5, label = 'Qhb')
a2.plot(time, sl3, 'b', linewidth = 1.5, label = 'Fbc')
a2.plot(time, sl4, 'g', linewidth = 1.5, label = 'S')
a2.set_ylabel('Traditional Sensitivity functions')
a2.legend(loc = 'best')
a2.set_xlabel('time (min)')

#compute generalise sensitivity
#2 by 2 parameter are taken into account
N1 = 2
V1 = sp.zeros((len(sl1),N1))
V1[:,0] = sl1
V1[:,1] = sl2

N2 = 2
V2 = sp.zeros((len(sl1),N2))
V2[:,0] = sl1
V2[:,1] = sl3

N3 = 2
V3 = sp.zeros((len(sl1),N3))
V3[:,0] = sl1
V3[:,1] = sl4

N4 = 2
V4 = sp.zeros((len(sl1),N4))
V4[:,0] = sl2
V4[:,1] = sl3

N5 = 2
V5 = sp.zeros((len(sl1),N5))
V5[:,0] = sl2
V5[:,1] = sl4

N6 = 2
V6 = sp.zeros((len(sl1),N6))
V6[:,0] = sl3
V6[:,1] = sl4


#matrix
D = sp.matrix(sp.zeros((N1,N1)) )
Dinv = sp.matrix(sp.zeros((N1,N1)) )

D2 = sp.matrix(sp.zeros((N2,N2)) )
Dinv2 = sp.matrix(sp.zeros((N2,N2)) )

D3 = sp.matrix(sp.zeros((N3,N3)) )
Dinv3 = sp.matrix(sp.zeros((N3,N3)) )

D4 = sp.matrix(sp.zeros((N4,N4)) )
Dinv4 = sp.matrix(sp.zeros((N4,N4)) )

D5 = sp.matrix(sp.zeros((N5,N5)) )
Dinv5 = sp.matrix(sp.zeros((N5,N5)) )

D6 = sp.matrix(sp.zeros((N6,N6)) )
Dinv6 = sp.matrix(sp.zeros((N6,N6)) )

for i in range(len(sl1)) :
    tmp1 = sp.matrix(V1[i,:])
    tmp2 = sp.matrix(V2[i,:])
    tmp3 = sp.matrix(V3[i,:])
    tmp4 = sp.matrix(V4[i,:])
    tmp5 = sp.matrix(V5[i,:])
    tmp6 = sp.matrix(V6[i,:])
    
    D = D + tmp1.T * tmp1
    D2 = D2 + tmp2.T * tmp2
    D3 = D3 + tmp3.T * tmp3

    D4 = D4 + tmp4.T * tmp4
    D5 = D5 + tmp5.T * tmp5
    D6 = D6 + tmp6.T * tmp6

Dinv = np.linalg.inv(D)
Dinv2 = np.linalg.inv(D2)
Dinv3 = np.linalg.inv(D3)
Dinv4 = np.linalg.inv(D4)
Dinv5 = np.linalg.inv(D5)
Dinv6 = np.linalg.inv(D6)

GS = sp.zeros((N1, len(sl1)))
G = sp.zeros((N1, 1))

GS2 = sp.zeros((N2, len(sl1)))
G2 = sp.zeros((N2,1))

GS3 = sp.zeros((N3, len(sl1)))
G3 = sp.zeros((N3,1))

GS4 = sp.zeros((N4, len(sl1)))
G4 = sp.zeros((N4, 1))

GS5 = sp.zeros((N5, len(sl1)))
G5 = sp.zeros((N5,1))

GS6 = sp.zeros((N6, len(sl1)))
G6 = sp.zeros((N6,1))

for i in range(len(sl1)) :
    Grad_L = sp.matrix(V1[i,:]).T
    Grad_L2 = sp.matrix(V2[i,:]).T
    Grad_L3 = sp.matrix(V3[i,:]).T
    Grad_L4 = sp.matrix(V4[i,:]).T
    Grad_L5 = sp.matrix(V5[i,:]).T
    Grad_L6 = sp.matrix(V6[i,:]).T

    tmp1 = Dinv * Grad_L
    tmp2 = Dinv2 * Grad_L2
    tmp3 = Dinv3 * Grad_L3
    tmp4 = Dinv4 * Grad_L4
    tmp5 = Dinv5 * Grad_L5
    tmp6 = Dinv6 * Grad_L6

    G = G + np.multiply(tmp1,Grad_L)
    G2 = G2 + np.multiply(tmp2,Grad_L2)
    G3 = G3 + np.multiply(tmp3,Grad_L3)
    G4 = G4 + np.multiply(tmp4,Grad_L4)
    G5 = G5 + np.multiply(tmp5,Grad_L5)
    G6 = G6 + np.multiply(tmp6,Grad_L6)

    GS[:,i] = G.T
    GS2[:,i] = G2.T
    GS3[:,i] = G3.T
    GS4[:,i] = G4.T
    GS5[:,i] = G5.T
    GS6[:,i] = G6.T

   
#Plot all :
Fig4 = plt.figure(figsize = [15.,10.])
a1 = Fig4.add_subplot(231)
a1.set_title('A')
a2 = Fig4.add_subplot(232)
a2.set_title('B')
a3 = Fig4.add_subplot(233)
a3.set_title('C')
a4 = Fig4.add_subplot(234)
a4.set_title('D')
a5 = Fig4.add_subplot(235)
a5.set_title('E')
a6 = Fig4.add_subplot(236)
a6.set_title('F')


a1.plot(time, GS[0,:], 'k', linewidth = 1.5, label = 'Ksh')
a1.plot(time, GS[1,:], 'r', linewidth = 1.5, label = 'Qhb')
a1.set_ylabel('Generalized Sensitivity functions')
a1.legend(loc = 'best')
a1.set_xlabel('time (min)')

a2.plot(time, GS2[0,:], 'k', linewidth = 1.5, label = 'Ksh')
a2.plot(time, GS2[1,:], 'b', linewidth = 1.5, label = 'Fbc')
a2.set_ylabel('Generalized Sensitivity functions')
a2.set_xlabel('time (min)')
a2.legend(loc = 'best')

a3.plot(time, GS3[0,:], 'k', linewidth = 1.5, label = 'Ksh')
a3.plot(time, GS3[1,:], 'g', linewidth = 1.5, label = 'S')
a3.set_xlabel('time (min)')
a3.set_ylabel('Generalized Sensitivity functions')
a3.legend(loc = 'best')

a4.plot(time, GS4[0,:], '-r', linewidth = 1.5, label = 'Qhb')
a4.plot(time, GS4[1,:], '-b', linewidth = 1.5, label = 'Fbc')
a4.set_ylabel('Generalized Sensitivity functions')
a4.legend(loc = 'best')
a4.set_xlabel('time (min)')

a5.plot(time, GS5[0,:], '-r', linewidth = 1.5, label = 'Qhb')
a5.plot(time, GS5[1,:], '-g', linewidth = 1.5, label = 'S')
a5.set_ylabel('Generalized Sensitivity functions')
a5.legend(loc = 'best')
a5.set_xlabel('time (min)')

a6.plot(time, GS6[0,:], '-b', linewidth = 1.5, label = 'Fbc')
a6.plot(time, GS6[1,:], '-g', linewidth = 1.5, label = 'S')
a6.set_ylabel('Generalized Sensitivity functions')
a6.legend(loc = 'best')
a6.set_xlabel('time (min)')


plt.show()
