#!/usr/bin/env python
# check prediction of P(k) for mocks
# generate spectra with pixel dr, apply Gaussian smearing
# and compute P(k) using bins that re smaller by a factor nSplit.
# when nSplit goes to infinity the measured P(k) should agree with the prediction
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt 
from scipy import signal

PI = np.pi

NX = 100
DX = 2. # Mpc/h
LX = NX * DX 
k_ny = PI /DX 
nkBin = NX//2 + 1   # NX even
nSpec = 100000
direct = True

Nsplit = 10
DXs =  DX / Nsplit
NXs = NX * Nsplit
k_nys = PI /DXs

smearing = True
nBinWin = 3 # in DX units: 3 in the mocks, NX if we want full Gaussian
nBinWin *= Nsplit
sigma = DX # in Mpc/h
sigma *= Nsplit / DX  # in units of DX / Nsplit

#********************************************************************
def P0(k) :     # find a function that is defined at k=0   <==
    k = np.abs(k)
    cut = np.not_equal(k,0) 
    eps = 1E-20
    k = k + eps * np.equal(k,0)
    P = -3.4 + 0.175 -0.35*np.log(k)-0.075*np.log(k)*np.log(k);
    P = PI * np.exp(P) /k;  
    #P[...]=1	# prov
    return P * cut 	# P(k)=0 for  k=0 

#********************************************************************
def Pcut(k,epsilon=0) :		# cut P(k) above k_ny
  P = P0(k)
  cut = np.less_equal(np.abs(k),k_ny-epsilon)
  return P * cut 	# P(k)=0 for  k > k_ny


#********************************************************************
def GenerateSpectrumDirect(N_X,weight):
    nkBin = N_X//2 + 1   # NX even
    xx = weight * sp.random.standard_normal(nkBin) 	
    yy = weight * sp.random.standard_normal(nkBin) 	
    z = (xx + 1j * yy)/np.sqrt(2)
    z[0] = xx[0]
    if N_X%2 == 0 : z[N_X//2] = xx[N_X//2]
    delta = np.fft.irfft(z) * np.sqrt(N_X)  # to get normalized iDFT
    return delta 

#********************************************************************
def SplitPixel(delta,k) :   # seems to be faster than np.repeat which does exactly the same !!   <==
    # split eack pixel by a factor k, giving the same original values to the resulting k pixels 
    # e.g. for k=4:  (x,y,z) -> (x,x,x,x,y,y,y,y,z,z,z,z)
    y = np.ones((k,len(delta)))
    z = y * delta
    delta2 = (z.T).reshape(k*len(delta))
    return delta2


#*************************************************************
def Prediction(q1,q2,k) :
    P = np.zeros(len(k))
    for q in np.arange(q1,q2) :
        W = sinc(k*DX/2)
        if (smearing): W*= np.exp(-k*k*DX*DX/2)
        if (q%2==0):
            P += Pcut(k+2*q*k_ny,epsilon=-0.000001)*W*W
        else:
            P += Pcut(k+2*q*k_ny,epsilon=0.000001)*W*W
    return P


#*************************************************************
#   average m pixels together, leaving away the N%m remaining pixels
def regroup(spectrum,m) :
    m = int(m)
    N=len(spectrum)
    p = N / m
    spectrum = spectrum[0:p*m]
    xx = spectrum.reshape(p,m)
    xx = xx.mean(1)
    return xx

#*************************************************************
def sinc(x) :
    return np.sinc(x/PI)  # np.sinx(x)=sin(PI*x)/(PI*x)

#********************************************************************
# 1D FFT for an even real function
# in:  a(0, ....,xmax)
# out: b[0,...,kmax] , 
def erFFT(a) :
  a_inv = a[-2:0:-1]    # a[-2], ... , a[1]
  a = np.append(a,a_inv)
  b_cmplx = np.fft.rfft(a)
  b = b_cmplx.real  	
  return  b





#*************************************************************
#*************************************************************      main
#.......................   compute weight(k)
k = np.fft.rfftfreq(NX) * 2 * k_ny	# rfftfreq 0 -> 0.5
P = P0(k)
weight = np.sqrt(P/DX)
ks = np.fft.rfftfreq(NXs) * 2 * k_nys	# rfftfreq 0 -> 0.5
Ps = P0(ks)

print "k_N=", k_ny

if (smearing) :
    window = signal.windows.gaussian(nBinWin,sigma)
Pmeas = np.zeros(len(ks))
for i in range(nSpec) :	#------------------------------ loop on spectra --------
    delta = GenerateSpectrumDirect(NX,weight) 
    delta = SplitPixel(delta,Nsplit)
    #print(len(delta))
    #delta  = np.repeat(delta,Nsplit)
    if (smearing) : 
        delta = np.concatenate((delta,delta[0:nBinWin-1]))
        delta = signal.convolve(delta, window, mode='valid') / sum(window)
    #print(len(delta))
    Gamma = np.fft.rfft(delta)
    #print(len(Gamma),len(Pmeas))
    Pmeas += (DXs/NXs)*abs(Gamma)**2  
Pmeas /= nSpec

print "ploting"
pred = Prediction(-Nsplit/2,Nsplit/2,ks)
plt.plot(ks/k_ny,Pmeas,color="blue")  
plt.plot(ks/k_ny,P0(ks),color="red")
plt.plot(ks/k_ny,pred,color="green")
plt.grid()
plt.yscale('log')
plt.show()

cut=np.where(pred>1E-30)
#print cut
plt.plot(ks[cut]/k_ny,(pred[cut]-Pmeas[cut])/pred[cut],color="green")
plt.grid()
plt.show()


xi0 = erFFT(Ps)
xiMeas = erFFT(Pmeas)
r=np.arange(len(xi0))* DX+NX/len(xi0)
plt.plot(r,r*xi0)
plt.plot(r,r*xiMeas)
plt.plot(r,r*(xiMeas-xi0))
plt.grid()
plt.show()
