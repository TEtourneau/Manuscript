'''  Fit P1D data over 2.2 < z < 3.6
produce a model that use Arinyo-i-Prats to extrapolate outside 0.2 < k < 2.0
extrapolation to z =2 might be done as P(2) = P(2.2) * P(2.2)/P(2.4) 
'''
from __future__ import division, print_function
from SaclayMocks import powerspectrum, util, constant
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt 
from astropy.cosmology import Planck15 as cosmo
from fitsio import FITS
import iminuit

kms = False     # km/s or Mpc/h
plotratio = False    #  plot P1D(k,z) / P1D(k,z0)
fitOld = False      # fitting P1D data, obsolete ?
kPk = False      # plot k P(k)/pi or just P(k)
fitSi = True
fitSiAll = False    # not working yet

def ReadPk(zref,kms=True):
    if (kms): xx = 1
    else: xx = cosmo.H(zref).value/(1+zref)/cosmo.h   # km/s -> Mpc/h
    msk = np.where( np.abs(z - zref) < 0.001)[0]
    if len(msk) == 0:
        # z from 2.2 to 4.4
        raise ValueError("ERROR -- You entered a wrong redshift. Here is the list of redshift : {}".format(np.unique(z)))
    k= data[:, 1][msk]*xx
    Pk= data[:, 2][msk]/xx
    Pkerr = data[:, 3][msk]/xx
    return k, Pk, Pkerr

#data = np.loadtxt("../etc/pk_fft35bins_noSi.out") # DR9
# data = np.loadtxt("../etc/pk_1d_DR12_13bins_noSi.out") # DR12 no Si III
data = np.loadtxt("/Users/tetourne/Work/Thesis/Git/SaclayMocks/etc/pk_1d_DR12_13bins_noSi.out") # DR12 no Si III
#data = np.loadtxt("../etc/pk_1d_DR12_13bins_no_BAL_DLA_newDLAmask.out") # DR12 no Si III
z = data[:, 0]
zmin=2.2
zmax=3.6
nz = int(round(1+(zmax-zmin)/0.2))

if (plotratio):  # plot P1D(k,z) / P1D(k,z0)
    from scipy.optimize import curve_fit
    def ff(z,a,gamma):
        return a * (1+z)**gamma
    yy = 1.2167
    _, Pk0, _ = ReadPk(2.8,kms)  # smallest errors at 2.8
    
    means = np.zeros(nz)
    errors = np.zeros(nz)
    zz = np.linspace(zmin,zmax,nz)
    print(zz)
    for i,zref in enumerate( zz) :
        #if i in [0,1,2,3,4,5,6,7]:  # which z do you want to plot
        if True:  # which z do you want to plot
            k, Pk, Pkerr = ReadPk(zref,kms)
            cor = 1 / yy**i #  to get all curves on top of each other
            cor /= Pk0 / yy**3 #   to get all curves around 1
            cor = 1 / Pk0   # to get all curves nearly flat
            Pk *=cor
            Pkerr *=cor
            #plt.errorbar(k*(1+0.003*i/k),cor*Pk,yerr=cor*Pkerr)
            plt.errorbar(k,Pk,yerr=Pkerr)
            weights = 1/Pkerr**2
            #weights /= weights/10
            means[i] = np.average(Pk,weights=weights)
            errors[i] = 1 / np.sqrt(weights.sum()) 
            print(i,means[i],errors[i])
            plt.plot(k, means[i]*np.ones(len(Pk)))
    plt.show()
    plt.errorbar(zz,means,yerr=errors,fmt="o")
    #plt.xscale("log")
    #plt.yscale("log")
    p, err = curve_fit(ff,zz,means,sigma=errors)
    print (p[0], np.sqrt(err[0,0]) )
    print (p[1], np.sqrt(err[1,1]) )
    zz=np.linspace(zz[0],zz[-1],100)
    Amp = ff(zz,p[0],p[1])
    plt.plot(zz,Amp)
    plt.show()

def hfitga(plt,xFit,yFit,error,color,label):
    
    a = -61.3
    b = 3
    c = 0
    
    def g(k,a,b,c):
        return np.exp(a*k+b)+c 
     
    def chi2TestM(a,b,c):    
        return np.sum((yFit-g(xFit,a,b,c))**2/error**2)
    
    fit_arg = dict( a =a, limit_a =[-100.,0.],   error_a =0.1,
                    b=b,limit_b=[0.,10],   error_b=0.04,
                    c =c, limit_c =[0.,100.], error_c =1.)
    m=iminuit.Minuit(chi2TestM,errordef=1.,**fit_arg)#,print_level=1
    m.migrad()
    print('******** migrad *******************************')
    print('      chi2=',m.fval)
    print('a=',m.values['a'],'+/-',m.errors['a'])
    print('b=',m.values['b'],'+/-',m.errors['b'])
    print('c=',m.values['c'],'+/-',m.errors['c'])
#    print '***********************************************'
#    plt.ylim(0.5,m.values['norm']*1.2)
    plt.errorbar(xFit,yFit,error,color="blue",label=label)
    x=np.linspace(xFit[0],xFit[-1],200)
    plt.plot(x,g(x,m.values['a'],m.values['b'],m.values['c']),color=color,linewidth=2)
#    plt.yscale('log')
    plt.show()
    return m

if (False):
    def g(k,a,b):
        return np.exp(a*k+b) 
    def chi2TestM(a,b):  
        return np.sum(( Pk-g(k,a,b) )**2/Pkerr**2)
    k, Pk, Pkerr = ReadPk(2.2,kms)
    #hfitga(plt,k,Pk,Pkerr,"red","essai")
    #exit(0)
    a=-61.3
    b=3
    fit_arg = dict( a = a, limit_a=[-100,0] ,error_a = 0.1,
                    b = b, limit_b=[0,10] ,error_b = 0.04)
    m=iminuit.Minuit(chi2TestM,errordef=1.,**fit_arg)#,print_level=1 
    m.migrad()
    a = m.values['a']
    b = m.values['b']
    print(a,b)
    plt.errorbar(k,Pk,Pkerr,fmt='o')
    plt.plot(k,g(k,a,b),color="red")
    plt.show()
    exit(0)

if (False):
    import iminuit

    def chi2test(a,b,c,d):  
        return np.sum(( Pk-g(k,a,b,c,d) )**2/Pkerr**2)
    def g(k,a,b,c,d):
        return np.exp(a*k+b+c/(k+d)) 
    zz=2.2
    k, Pk, Pkerr = ReadPk(zz,kms)
    fit_arg = dict( a = -61.3, limit_a=[-300,0] ,error_a = 0.04,
                    b = 3, limit_b=[-10,10] ,error_b = 0.04,
                    c = 0, limit_c=[0,10] ,error_c = 0.04,
                    d = 0, limit_d=[-1,1] ,error_d = 0.04)
    m=iminuit.Minuit(chi2test,errordef=1.,print_level=3,**fit_arg)#,print_level=1 
    m.migrad()
    a = m.values['a']
    b = m.values['b']
    c = m.values['c']
    d = m.values['d']
    print(a,b,c,d)
    k, Pk, Pkerr = ReadPk(zz,kms)
    plt.errorbar(k,Pk,Pkerr,fmt='o')
    plt.plot(k,g(k,a,b,c,d))
    plt.show()
    exit(0)


if (fitSi):    # fit one z at a time, including fixed SiII and SiIII contributions from NPD
                # this fit must be done in km/s, once can then transform the parameters to Mpc/h
    data2=np.loadtxt("P1D_Prats.txt")
    k=data2[:,0]
    P=data2[:,1]
    Pprats = sp.interpolate.InterpolatedUnivariateSpline(k,P)
    import iminuit
    plotratio = False # plot z * data / fit 
    exp_cst = False # exp + cst : exp(a*k+b)+c
    nz = int(round((zmax-zmin)/0.2)) + 1
    teff = np.array([0.184924702397, 0.231425921518, 0.285929695332, 0.349252333224, 0.422242531447,0.505780851386, 0.600779232412, 0.708180535481, 0.828958114192, 0.9641154105, 1.1146855727, 1.28173109358,1.46634346696 ])
    zt = np.linspace(2.2,4.6,13)
    teffz =sp.interpolate.InterpolatedUnivariateSpline(zt,teff)
    def chi2test(a,b,c,d):  
        return np.sum(( Pk-g(k,a,b,c,d) )**2/Pkerr**2)
    def chi2f(a,b,c,d):  
        return np.sum(( Pk-f(k,a,b,c) )**2/Pkerr**2)
    def fSi(k):
        dv2 = 5577. # SiII
        f2 = 7.25404e-04
        A2 = f2/(1-np.exp(-teffz(zz)))
        dv3 = 2271. # SiIII
        f3 = 5.98278e-03
        A3 = f3/(1-np.exp(-teffz(zz)))
        #return  (1+A3*A3+A3*np.cos(k*dv3))
        return (1+A2*A2+A2*np.cos(k*dv2)) * (1+A3*A3+A3*np.cos(k*dv3))
    def g0(k,a,b,c,d):
        return np.exp(-a*k+b+c/(k+d))
    def g(k,a,b,c,d):
        return  g0(k,a,b,c,d) *  fSi(k)
    def f0(k,a,b,c):
        return np.exp(-a*k+b)+c
    def f(k,a,b,c):
        return  f0(k,a,b,c) *  fSi(k)
    
    outfits = FITS("P1Dmodel.fits", 'rw', clobber=True)
    #zarray = np.zeros(0); karray = np.zeros(0); Parray = np.zeros(0); 
    zarray = []; karray = []; Parray = []; 
    a=np.zeros(nz); b=np.zeros(nz); c=np.zeros(nz); d=np.zeros(nz); sig2LSS=np.zeros(nz)
    for i, zz in enumerate( np.linspace(zmin,zmax,nz) ):
        k, Pk, Pkerr = ReadPk(zz,True)  # km/s
        if (exp_cst) :  # exp(a*k+b)+c, does not go to zero at high k
            fit_arg = dict( a = 61.3, limit_a=[0,300] ,error_a = 0.04,
                        b = 3, limit_b=[-10,10] ,error_b = 0.04,
                        c = 0, limit_c=[0,10] ,error_c = 0.04)
            m=iminuit.Minuit(chi2f,errordef=1.,print_level=0,**fit_arg)#,print_level=1 
        else :      # exp(-a*k+b+c/(k+d))
            fit_arg = dict( a = 48,  limit_a=[0,300] ,error_a = 0.04, fix_a=True,
                        b = 3, limit_b=[-10,10] ,error_b = 0.04,
                        c = 0, limit_c=[0,0.005] ,error_c = 0.001,
                        d = 0, limit_d=[0,0.007] ,error_d = 0.001)
            m=iminuit.Minuit(chi2test,errordef=1.,print_level=0,**fit_arg)#,print_level=1 
        m.migrad()
        kms2Mpc = cosmo.H(zz).value/(1+zz)/cosmo.h   # ~ 100, to get param in Mpc/h
        #kms2Mpc = 1 # prov !!!!
        a[i] = m.values['a']/kms2Mpc
        b[i] = m.values['b']-np.log(kms2Mpc)
        c[i] = m.values['c']*kms2Mpc
        d[i] = m.values['d']*kms2Mpc
        #print(round(10*zz)/10,m.values['a'],m.values['b'],m.values['c'],m.values['d']) # in km/s for simu1D
        print(round(10*zz)/10,a[i],b[i],c[i],d[i]) # in Mpc/h
        #continue  # prov !!!
        #print(zz,a[i]/xx,b[i]-np.log(xx),xx*c[i],xx*d[i])
        #k, Pk, Pkerr = ReadPk(zz,kms)
        if (kPk) : yy = k/np.pi    # plot k P(k)/pi
        else: yy=1            # plot P(k)
        if (plotratio): yy = zz/g(k,a[i],b[i],c[i],d[i])
        k, Pk, Pkerr = ReadPk(zz,False) # Mpc/h
        plt.errorbar(k,yy*Pk,yy*Pkerr,fmt='o', color=constant.colors[i], label="z = {}".format(np.round(zz,1)))
        kk=np.linspace(0.0,20,1000)
        if (exp_cst): PP = yy*f(k,a[i],b[i],c[i])
        #else: PP = yy*g(k,a[i],b[i],c[i],d[i])
        else:  PP = yy*g0(kk,a[i],b[i],c[i],d[i])
        cut = (kk < 0.2)
        ii = len(kk[cut])
        PP[cut] = (PP[ii]/Pprats(kk[ii])) * Pprats(kk[cut]) # so PP[ii] is unchanged, continuous
        plt.plot(kk,PP, color=constant.colors[i])
        #if (i==2): #prov
        #    plt.plot(kk,PP)
        karray.append(kk)
        Parray.append(PP)
        zarray.append(zz*np.ones(len(kk))) 
        #........... compute sigma_F
        #kk /= kms2Mpc
        #PP *= kms2Mpc
        Pofk = sp.interpolate.InterpolatedUnivariateSpline(kk,PP)
        #sigLSS = util.sigma_p1d(Pofk,69)
        sigLSS = util.sigma_p1d(p1dmiss=Pofk)
        print (sigLSS)
        sig2LSS[i] = sigLSS**2
        #prov print(zz,a[i],b[i],c[i],d[i],sigLSS**2)
        if (False):  # check sigma_LSS
            kny = np.pi / 69
            cut = (kk<kny)
            genData = powerspectrum.GenerateData(kk[cut],PP[cut])
            deltas=[]
            for iii in np.arange(10000):
                delta, DX  = genData.oneArray()
                deltas.append(delta)
            deltas=np.concatenate(deltas)
            #print(DX,deltas.std()**2,sigLSS**2)
        if (exp_cst): plt.plot(k,yy*f(k,a[i],b[i],c[i]))
        #else: plt.plot(k,yy*g(k,a[i],b[i],c[i],d[i]))
        #else: plt.plot(kk,yy*g0(kk,a[i],b[i],c[i],d[i]))
    if (kPk): plt.yscale("log")
    karray=np.concatenate(karray)
    Parray=np.concatenate(Parray)
    zarray=np.concatenate(zarray)
    outfits.write(zarray, extname="z")
    outfits.write(karray, extname="k")
    outfits.write(Parray, extname="P")
    outfits.close()
    plt.yscale("log")
    #for i in range(nz): print(zmin+0.2*i,a[i],b[i],c[i],d[i])
    plt.show()
    z = np.linspace(zmin,zmax,nz)
    l = 1216 * (1+z)
    #plt.plot(l,sig2LSS/2)
    # f=FITS("../essai/delta_attributes.fits")
    # data=f[2].read()
    # l=10**(data["loglam"])
    # varlss = data["var_lss"]
    #plt.plot(l,varlss)    
    #plt.show()

if (fitOld):  # fitting P1D data, obsolete ?
    from scipy.optimize import curve_fit
    
    def f(k,a,b):
        return np.exp(a*k+b)
    def f1(k,b):
        return np.exp(-67*k+b)    # km/s
        #return np.exp(-0.67*k+b)    # mpc/h
    def g(k,a,b,c,d):
        return np.exp(a*k+b+c/(k+d))
    #def g(k,a,b,c,d):
        #return np.exp(a*k+b+c/(k+d))
        #return np.exp(a*k+b+c*np.exp(-d*k))
    def g2(k,c,d):
        a=-0.55857373  # fit to Prats' model
        b=-1.21583873
        return np.exp(a*k+b+c/(k+d))
    
    for i,zref in enumerate( np.linspace(zmin,zmax,8)) :
        if i in [0,1,2,3,4,5,6,7]:  # which z do you want to plot
        #if i in [3]:  # which z do you want to plot
            if (kPk) : yy = k/np.pi    # plot k P(k)/pi
            else: yy=1            # plot P(k)
            k, Pk, Pkerr = ReadPk(zref,kms)
            if (kms): p0f= [-70,3]
            else : p0f= [-1,-2]; p0g=[-1,-2,0,0]
            #par, err =curve_fit(f1,k,Pk,p0=[-1],sigma=Pkerr)
            #plt.plot(k,yy*f1(k,par[0]))
            par, err =curve_fit(f,k,Pk,p0=p0f,sigma=Pkerr)
            plt.plot(k,yy*f(k,par[0],par[1]))
            #p0g[0:2]=par
            #print (p0g)
            #par, err =curve_fit(g,k,Pk,p0=p0g,sigma=Pkerr)
            #par, err =curve_fit(g2,k,Pk,p0=[0,0],sigma=Pkerr)
            print (zref,par)
            plt.errorbar(k,yy*Pk,yy*Pkerr,fmt='o')
            kk=np.linspace(0.1,15,1000)
            #plt.plot(kk,yy*g2(kk,par[0],par[1]),color="red")
            #plt.plot(kk,yy*g2(kk,0,0),color="blue")
            #plt.plot(kk,yy*f(kk,par[0],par[1],par[2]))
            plt.yscale("log")
    if (kPk): plt.yscale("log")
    plt.show()

'''

if (fitSiAll):    # fit all z simultaneously including SiII and SiIII contributions
    #   not working yet
    import iminuit
    zmin=2.2
    zmax=2.2
    nz = int(round((zmax-zmin)/0.2)) + 1
    teff = np.array([0.184924702397, 0.231425921518, 0.285929695332, 0.349252333224, 0.422242531447,0.505780851386, 0.600779232412, 0.708180535481, 0.828958114192, 0.9641154105, 1.1146855727, 1.28173109358,1.46634346696 ])
    zt = np.linspace(2.2,4.6,13)
    teffz =sp.interpolate.InterpolatedUnivariateSpline(zt,teff)
    def chi2(a,b,c,d,f2,f3):  
        chi_2=0
        for zz in np.linspace(zmin,zmax,nz):
            k, Pk, Pkerr = ReadPk(zz,kms)
            bb = b + 1.2 * (zz-zmin) /1.4
            chi_2 += ((Pk-g(k,a,bb,c,d,f2,f3))**2/Pkerr**2)
        return chi_2
    def g(k,a,b,c,d,f2,f3):
        dv3 = 2271. # SiIII
        dv2 = 5577. # SiII
        A2 = f2/(1-np.exp(-teffz(zz)))
        A3 = f3/(1-np.exp(-teffz(zz)))
        return np.exp(a*k+b+c/(k+d)) * (1+A2*A2+A2*np.cos(k*dv2)) * (1+A3*A3+A3*np.cos(k*dv3))
    
    fit_arg = dict( a = -61, error_a = 0.04,
                    b = 3, error_b = 0.04,
                    c = 0, error_c = 0.04,
                    d = 0, error_d = 0.04,
                    f2 = 0, error_f2 = 0.04,
                    f3 = 0, error_f3 = 0.04)
    # mean =mean, limit_mean =[-2.,2.],   error_mean =0.04,
    m=iminuit.Minuit(chi2,errordef=1.,**fit_arg)#,print_level=1 
    a = m.values['a']
    b = m.values['b']
    c = m.values['c']
    d = m.values['d']
    f2 = m.values['f2']
    f3 = m.values['f3']
    print(a,b,c,d,f2,f3)

    if (kPk) : yy = k/np.pi    # plot k P(k)/pi
    else: yy=1            # plot P(k)
    for zz in np.linspace(zmin,zmax,nz):
        k, Pk, Pkerr = ReadPk(zz,kms)
        plt.errorbar(k,yy*Pk,yy*Pkerr,fmt='o')
        plt.plot(k,g(k,a,b,c,d,f2,f3))
    plt.show()

'''
