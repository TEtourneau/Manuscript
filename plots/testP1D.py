#!/usr/bin/env python
# compare the prediction of P1D with/without sinc and Gaussian smearing
#
# should add: 
# test prediction of mock P1D(k)
# generate mocks with ~2 Mpc/h voxels, compute their P(k) using small pixels
# option to add Gaussian smearing
# compute corresponding CF, compare P(k) and CF to prediction
# with and without RSD
from __future__ import division, print_function
from SaclayMocks import powerspectrum, util
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt 
from scipy import interpolate

RSD = False # not working currently, see MissingPk.py   <==

PI = np.pi
dx = 2.19   # box voxel size
sigma =  1 * dx 	# in Mpc/h, if sigma <0 take a constant over large pixels
pixel = 0.2  # spectrum pixel
kmax = 20    # 20
dk=0.001    # with RSD option requires too much memory for laptop
dk=0.005;  # prov
k_ny = PI/dx

#********************************************************************
def sinc(x) :
    return np.sinc(x/PI)  # np.sinc(x)=sin(PI*x)/(PI*x)

#********************************************************************   to be put in powerspectrum.py <==
def xi_from_pk_1D(k,P1): # k must go from 0 to kmax with constant steps   to be improved <==
    kmax=np.max(k)
    if (kmax<100):  # zero padding
        #k
        PP=np.zeros(1)
    N = len(k)
    P2 = np.flipud(P1[1:-1])
    P=np.hstack((P1,P2))        #  P0,P1, .. PN,PN-1, ... P1
    xi=np.real( np.fft.rfft(P) )
    xi *= kmax/N
    rmax = np.pi * N / 2 / kmax
    r=np.linspace(0,rmax,N)
    return r,xi

#******************************************** Main() ************************

print ("k_ny=", k_ny, "sigma =", sigma)

#.................................................... define k vectors
kk=np.arange(kmax/dk)*dk   # kk to compute all P(k)
cut = (kk<=k_ny)
kk_cut=kk[cut]
if (RSD) :
    k_par = np.arange(kmax/dk)*dk
    k_par_t = k_par.reshape(len(k_par),-1)
    k_perp = np.arange(kmax/dk)*dk
    k_perp_t = k_perp.reshape(len(k_perp),-1)

kpmax = PI / pixel
nkbin = int(kpmax/ dk) + 1
kp = np.linspace(0,kpmax,nkbin)
print(nkbin,kp)

#kp = np.arange(PI/0.02/dk)*dk  # k for plot
#print(kp)


#.................................................... define P camb
P_camb = powerspectrum.P_0("$SACLAYMOCKS_BASE/etc/PlanckDR12.fits")
Pcamb = P_camb.P(kk)
Pcamb_cut=Pcamb[cut]

P1Dcamb = powerspectrum.P_1D(kk,Pcamb).P1D(kp)
plt.plot(kp,P1Dcamb,color="black")

#.................................     compute P1D with w^2 
# once we have W^2, the effect of k_ny cut is negligible
# so, indistinguishable from cut at k_N and W^2 
W = np.exp(- dx*dx*kk*kk/2)
P1DWcamb = powerspectrum.P_1D(kk,Pcamb*W*W).P1D(kp)
plt.plot(kp,P1DWcamb,color="blue")

if (RSD) :
    PRSD = P_RSD(k_par_t,k_perp,PW2)
    P1DWcambRSD = powerspectrum.P_1D_RSD(k_par,k_perp,PRSD).P1D(kp)
    plt.plot(kp,P1DWcambRSD,color="black")
#print ("1")

#.................................     compute P1D with cut at k_Nyquist
P1Dcutcamb = powerspectrum.P_1D(kk_cut,Pcamb_cut).P1D(kp)
#plt.plot(kp,P1Dcutcamb,color="blue")
plt.plot(kp,P1Dcutcamb,color="green")

if (RSD) :
    PRSD = P_RSD(k_par_t,k_perp,Pcut)
    P1DcutcambRSD = powerspectrum.P_1D_RSD(k_par,k_perp,PRSD).P1D(kp)
    plt.plot(kp,P1DcutcambRSD,color="green")
#print ("2")

#.................................     compute P1D with cut at k_N and W^2
W = np.exp(- dx*dx*kk_cut*kk_cut/2)
P1DWcutcamb = powerspectrum.P_1D(kk_cut,Pcamb_cut*W*W).P1D(kp)
plt.plot(kp,P1DWcutcamb,color="red")

if (RSD) :
    PRSD = P_RSD(k_par_t,k_perp,PW2cut)
    P1DWcutcambRSD = powerspectrum.P_1D_RSD(k_par,k_perp,PRSD).P1D(kp)
    plt.plot(kp,P1DWcutcambRSD,color="green")
#print ("3")

#.................................      missing P^1D(k)
P1Dmissing = np.maximum(interpolate.InterpolatedUnivariateSpline(kp,P1Dcamb - P1DWcutcamb),0)  
#plt.plot(kp,P1Dmissing(kp),color="red")

if (RSD) :
    P1DmissingRSD = np.maximum(interpolate.InterpolatedUnivariateSpline(kp,P1DcambRSD - P1DWcutcambRSD),0) 
    plt.plot(kp,P1DmissingRSD(kp),color="red")

#.................................      sigma's
if (kp.max() > 0. ):
    sig = util.sigma_p1d(p1dmiss=sp.interpolate.interp1d(kp,P1Dcamb),pixel=pixel)
    sigcut = util.sigma_p1d(p1dmiss=sp.interpolate.interp1d(kp,P1Dcutcamb),pixel=pixel)
    sigW = util.sigma_p1d(p1dmiss=sp.interpolate.interp1d(kp,P1DWcamb),pixel=pixel)
    sigWcut = util.sigma_p1d(p1dmiss=sp.interpolate.interp1d(kp,P1DWcutcamb),pixel=pixel)
    print("sig=",sig)
    print("sigcut=",sigcut)
    print("sigW=",sigW)
    print("sigWcut=",sigWcut)

print("plotting")
plt.grid()
plt.yscale("log")
plt.show()

#..................................... effect of cut at k_N after W_gaussian 
cut = np.where(P1DWcamb > 0)
plt.plot(kp[cut],P1DWcutcamb[cut]/P1DWcamb[cut])
plt.grid()
plt.show()


exit(0)

Delta_k = 0.0001

nkbin = int(k_ny / Delta_k) + 1
k1 = np.linspace(0,k_ny,nkbin)
k1[0] = 1E-20
Pk0 = P_camb.P(k1)
plt.plot(k1,Pk0,color='b')  # Camb

nkbin = int(kmax / Delta_k) + 1
k2 = np.linspace(0,kmax,nkbin)
k2[0] = 1E-20
Pk = P_camb.P(k2)
#plt.plot(k,Pk,'--',color='y')
PkW = P_camb.P(k2)*(np.exp(-k2*k2*sigma*sigma/2.))**2
if (sincfactor >0) : PkW *= sinc(k2*sincfactor*sigma)**2
Pk = P_camb.P(k2)*(np.exp(-k2*k2*sigma*sigma/2.))**2 * sinc(k2*sigma)**2 # prov
plt.plot(k2,Pk,color='r') # prov
plt.plot(k2,PkW,color='g')
plt.yscale('log')
plt.show()

#................................................
plt.xlim([0.,500 ])     # prov
plt.ylim([-30.,150 ])     # prov
#print (k2.max())
#r, xi0 = powerspectrum.xi_from_pk(k2,Pk)
r, xi0 = xi_from_pk_1D(k2,Pk)
print ("xi0 done")
plt.plot(r,r*r*xi0)
#r, xi = powerspectrum.xi_from_pk(k2,PkW)
r, xi = powerspectrum.xi_from_pk_1D(k2,PkW)
plt.plot(r, r*r*xi)
plt.plot(r, r*r*(xi-xi0))
plt.plot(r,0*r,color="black",ls="--",lw=0.5)
print ("ploting xi")
plt.show()

