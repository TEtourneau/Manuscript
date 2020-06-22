#!/usr/bin/env python
# effect of bining and Gaussian smearing on 3D P(k) and xi(r)   
# the effect of binning is approximated as isotropic
from __future__ import division, print_function
from SaclayMocks import powerspectrum, util
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt 
#from astropy.cosmology import Planck15 as cosmo

compareWsinc = True
compare_sincfactor = False
save = False

PI = np.pi
dx = 2.19
nPixel = 1024
LX = nPixel * dx
sigma =  1 * dx 	# in Mpc/h, if sigma <0 take a constant over large pixels
# sincfactor = np.sqrt(3) # if <0 no sinc factor srt(3)=1.732 is an overestimate
sincfactor = 1 # is an underestimate
# sincfactor = pow(6/np.pi,1./3)  # = 1.241, such that 4 pi r^3 / 3 = (2a)^3
kmax = 10.
Delta_k = 0.0001
k_ny = PI/dx

#********************************************************************
def sinc(x) :
    return np.sinc(x/PI)  # np.sinc(x)=sin(PI*x)/(PI*x)


#*********************************************   main   *********************

print ("k_ny=", k_ny, "sigma =", sigma)

#................ P_Camb
P_camb = powerspectrum.P_0("$SACLAYMOCKS_BASE/etc/PlanckDR12.fits")
nkbin = int(kmax / Delta_k) + 1
kfull = np.linspace(0,kmax,nkbin)
kfull[0] = 1E-20
Pkcambfull = P_camb.P(kfull)
# plt.plot(kfull,Pkcambfull,color='black',ls="--")  # Camb full
nkbin = int(k_ny / Delta_k) + 1
kcut = np.linspace(0,k_ny,nkbin)
kcut[0] = 1E-20
Pkcambcut = P_camb.P(kcut)
# plt.plot(kcut,Pkcambcut,color='black',)  # CambCut

#.................
WGaussian = (np.exp(-kfull*kfull*sigma*sigma/2.))**2
Wsinc = sinc(kfull*sincfactor*sigma/2)**2
PkW = Pkcambfull * WGaussian
Pksinc = Pkcambfull * Wsinc
PkWsinc = Pkcambfull * WGaussian * Wsinc
if (compareWsinc):
    # plt.plot(kfull,PkW,color="green")
    # plt.plot(kfull,PkWsinc,color="red")
    # plt.yscale("log")
    # plt.show()
    
#................................................ compute and plot xi
    plt.xlim([0.,200 ])     # prov
    plt.ylim([-20.,60 ])     # prov
    #print (kfull.max())
    r, xi0 = powerspectrum.xi_from_pk(kfull,Pkcambfull)
    print ("xi0 done")
    plt.plot(r,r*r*xi0,color="black")
    r, xi1 = powerspectrum.xi_from_pk(kfull,PkW)
    plt.plot(r, r*r*xi1,color="green")
    r, xi2 = powerspectrum.xi_from_pk(kfull,PkWsinc)
    plt.plot(r, r*r*xi2,color="red")
    
    plt.plot(r, r*r*(xi1-xi0),color="green")
    plt.plot(r, r*r*(xi2-xi0),color="red")
    plt.plot(r, r*r*(xi2-xi1),color="blue")
    plt.plot(r, 0*r,color="black",ls="--",lw=0.5)
    print ("ploting xi")
    plt.grid()
    plt.show()

#................................................ compare different sincfactor

if (compare_sincfactor) :
    sincfactor = 1
    Wsinc = sinc(kfull*sincfactor*sigma)**2
    PkWsinc1 = Pkcambfull * WGaussian * Wsinc
    sincfactor = np.sqrt(3)
    Wsinc = sinc(kfull*sincfactor*sigma)**2
    PkWsinc3 = Pkcambfull * WGaussian * Wsinc
    plt.plot(kfull,PkWsinc1,color="green")
    plt.plot(kfull,PkWsinc,color="blue")
    plt.plot(kfull,PkWsinc3,color="red")
    plt.yscale("log")
    plt.show()
    
    plt.xlim([0.,200 ])     # prov
    plt.ylim([-20.,60 ])     # prov
    r, xi1 = powerspectrum.xi_from_pk(kfull,PkWsinc1)
    plt.plot(r,r*r*xi1,color="green")
    r, xi = powerspectrum.xi_from_pk(kfull,PkWsinc)
    plt.plot(r, r*r*xi,color="blue")
    r, xi3 = powerspectrum.xi_from_pk(kfull,PkWsinc3)
    plt.plot(r, r*r*xi3,color="red")
    #plt.show()

    plt.plot(r, r*r*(xi1-xi),color="green")
    plt.plot(r, r*r*(xi3-xi),color="red")
    plt.plot(r, 0*r,color="black",ls="--",lw=0.5)
    print ("ploting xi")
    plt.show()
   


#................................................   voir save/test3D.py
if (save) :
    np.savez("mockPred3D.npz",k1=kcut,Pk0=Pkcambcut,kfull=kfull,PkW=PkW,r=r,xi=xi,xi0=xi0)
    npzfile=np.load("mockPred3D.npz")
    k1 = npzfile['k1']
    kfull = npzfile['kfull']
    Pk0 = npzfile['Pk0']
    PkW = npzfile['PkW']
    xi0 = npzfile['xi0']
    xi = npzfile['xi']
    r = npzfile['r']
    
