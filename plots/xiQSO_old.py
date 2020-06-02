'''
Coles and Jones (1991) show that if g Gaussian(mu,sigma) with xi_g
then exp(g) has CF = exp (xi_g) - 1
In order to have a z dependent QSO CF 
    we can draw QSO proportionally to exp[a(z) g] 
    with a(z_0)=1 and xi_QSO(z_0) = b_0^2 xi_m(z_0) = xi_0
    a^2 = b/(1+z) / [ b_0/(1+z_0) ] 
so xi_QSO = exp[ a^2 xi ] - 1  and exp [ xi ] - 1 = xi_0
    =>  xi_QSO = exp [ a^2 ln(1+xi_0)] - 1
if xi_0 << 1 or |a-1| << 1 then xi_QSO = a^2 xi_0 
Here we check how close we are. 

Then we go one step further, consider that we have two lognormal fields 
at z1 and z2. For a given z, we compute the probability with z_0 = z_1 and z_2
and interpolate linearly in z.

If we interpolate we are below, while we are above if we extrapolate.
So we can go even one step further and use 3 lognormal fields
ans make an average of interpolation and extrapolation to cancel the effect at r=5
The effect at other r is very small.

'''

import numpy as np
from SaclayMocks import powerspectrum
import matplotlib.pyplot as plt 

compareBiases = False
oneField = False
twoField = True
threeField = False
printing = False   # print to screen, should go to a file  <==

zref=2.33
growth = 0.37   # should compute growth(zref) instead <==
p=1     # plot r^p xi

def bias(z):
    #return 0.278*(1+z)**2 + 0.57
    return 3.7 * ((1+z)/(1+2.33))**1.7

def delta(z,z0):   # delta(z)/delta(z_0)
    return bias(z)/(1+z)/ (bias(z0)/(1+z0))

def plotxi_a(z,z0,r,xi0):
    xi1 = xi_a(z,z0,r,xi0)
    plt.plot(r,xi1*r**p)

def xi_a(z,z0,r,xi0):
    a=delta(z,z0)
    xi1 = np.exp( a*a *np.log(1+xi0))-1
    return xi1 / (a*a)

def xi_interp(z,z1,z2,r,xi0):
    xi1 = xi_a(z,z1,r,xi0)
    xi2 = xi_a(z,z2,r,xi0)
    return (xi1*(z2-z)+xi2*(z-z1))/(z2-z1)




if (compareBiases) :
    z=np.linspace(1.3,3.6,100)
    biasPL = bias(z)
    bH = 1.24 * ((1+z)/(1.59))**1.44
    bb = 3.7 * ((1+z)/3.33)**1.7
    plt.plot(z,bb/bH)
    plt.plot(z,bH)
    plt.plot(z,bb)
    plt.grid()
    plt.show()
    exit(0)


#................ compute xi_0 = xi_ln (zref)
kmax = 100
kk=np.linspace(0,kmax,100000)
# P_camb_ln = powerspectrum.P_ln("/Users/legoff/ana/LyaMocks/bin/data/PlanckDR12.fits",G_times_bias=growth*bias(zref))
P_camb_ln = powerspectrum.P_ln("~/Work/Thesis/Git/SaclayMocks/etc/PlanckDR12.fits",G_times_bias=growth*bias(zref))
Pln = P_camb_ln.P(kk)
r,xi0 = powerspectrum.xi_from_pk(kk,Pln)
xi0=xi0[r<200]
r=r[r<200]

if (oneField) :
    #z=np.linspace(2,5,100)
    #a=delta(z,z0)
    #plt.plot(z,a)
    #plt.show()
    z0 = 2.33
    xi2 = xi_a(2.0,z0,r,xi0)
    xi233 = xi_a(2.33,z0,r,xi0)
    xi26 = xi_a(2.6,z0,r,xi0)
    xi36 = xi_a(3.6,z0,r,xi0)
    xi455 = xi_a(4.55,z0,r,xi0) # z=4.55 => a=1.5
    plt.plot(r,xi2*r**p, label=r'$z=2$')
    plt.plot(r,xi26*r**p, label=r'$z=2.6$')
    plt.plot(r,xi36*r**p, label=r'$z=3.6$')
    plt.plot(r,xi233*r**p, label=r'$z=z_0$')
    plt.grid()
    plt.xlabel(r'$r$ [Mpc/h]')
    plt.ylabel(r'$r \xi(r)$')
    plt.legend()

    plt.figure()
    plt.plot(r,xi2/xi233, label=r'i : $z=2$')
    plt.plot(r,xi26/xi233, label=r'i : $z=2.6$')
    plt.plot(r,xi36/xi233, label=r'i : $z=3.6$')
    plt.grid()
    plt.xlabel(r'$r$ [Mpc/h]')
    plt.ylabel(r'$\xi_i / \xi_0$')
    plt.show()
    exit(0)

if (twoField) :
    z1 = 2.1
    z2 = 3.5
    for ztest in [1.9,2.8,3.6] :
#     for ztest in [2.1, 3.5] :
        xi=xi_interp(ztest,z1,z2,r,xi0)
        plt.plot(r,xi*r**p, label='z = {}'.format(ztest))
    plt.plot(r,xi0*r**p, label='z = 2.33')
    plt.grid()
    plt.xlabel(r'$r$ [Mpc/h]')
    plt.ylabel(r'$\xi_i / \xi_0$')
    plt.legend()
    exit(0)

if (False) :

    z1=1.9
    z2=2.75
    z3=3.6
    r0=5
    i0=len(r[r<r0])     # r(i0) ~ r0

    dz=0.01
    for ztest in np.linspace(z1+dz,z3-dz,int((z3-z1)/dz-1)) :  # to compute a at each z
    #for ztest in np.linspace(2.3,2.4,1) :   #  compute only for ztest = 2.3 to make a plot
        if (ztest==z2) : continue
        xi1=xi_interp(ztest,z1,z2,r,xi0)
        xi2=xi_interp(ztest,z2,z3,r,xi0)
        x1=xi1[i0]/xi0[i0]
        x2=xi2[i0]/xi0[i0]
        a = (1-x2)/(x1-x2)  # we want a x1 + (1-a) x2 = 1
        xi3=a*xi1+(1-a)*xi2
        if (ztest < z2) :
            print(ztest,a,max(xi3/xi0),min(xi3[i0:-1]/xi0[i0:-1]))
        else:
            print(ztest,a,min(xi3/xi0),max(xi3[i0:-1]/xi0[i0:-1]))
    
    plt.plot(r,xi1/xi0)
    plt.plot(r,xi2/xi0)
    plt.plot(r,xi3/xi0)
    plt.grid()
    plt.show()

if (threeField) :

    z1=1.9
    z2=2.75
    z3=3.6
    r0=5
    i0=len(r[r<r0])     # r(i0) ~ r0

    dz=0.01
    z = np.linspace(z1+dz,z3-dz,int((z3-z1)/dz-1))
    z = z[np.where(z!=z2)]
    a = np.zeros(len(z))
    for i,ztest in enumerate(z) :
        if (ztest==z2) : continue
        xi1=xi_interp(ztest,z1,z2,r,xi0)
        xi2=xi_interp(ztest,z2,z3,r,xi0)
        x1=xi1[i0]/xi0[i0]
        x2=xi2[i0]/xi0[i0]
        a[i] = (1-x2)/(x1-x2)  # we want a x1 + (1-a) x2 = 1
        xi3=a[i]*xi1+(1-a[i])*xi2
        if (ztest==2.3) :
            plt.plot(r,xi1/xi0)
            plt.plot(r,xi2/xi0)
            plt.plot(r,xi3/xi0)
        if (printing) :
            if (ztest < z2) :
                print(ztest,a[i],max(xi3/xi0),min(xi3[i0:-1]/xi0[i0:-1]))
            else:
                print(ztest,a[i],min(xi3/xi0),max(xi3[i0:-1]/xi0[i0:-1]))
    plt.grid()
    plt.show()
    plt.plot(z,a)
    plt.grid()
    plt.show()

exit(0)

#.................................   the use of code below is not too clear
xi1=xi_interp(2.85,z1,z2,r,xi0)
plt.plot(r,xi1*r**p)
plt.plot(r,xi0*r**p,color="orange")
xi2=xi_interp(2.8,z1,z2,r,xi0)
plt.plot(r,xi2*r**p,color="green")
xi3=xi_interp(2.9,z1,z2,r,xi0)
plt.plot(r,xi3*r**p,color="red")
xi4=xi_interp(2.7,z1,z2,r,xi0)
plt.plot(r,xi4*r**p,color="purple")
plt.grid()
plt.show()
plt.plot(r,xi1/xi0,color="orange")
plt.plot(r,xi2/xi0,color="green")
plt.plot(r,xi3/xi0,color="red")
plt.plot(r,xi4/xi0,color="purple")
plt.grid()
plt.show()
