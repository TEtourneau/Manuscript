import numpy as np
from SaclayMocks import powerspectrum
import matplotlib.pyplot as plt 


p = 2
d_cell = 2.19
kmax = 10
Delta_k = 0.0001
nkbin = int(kmax / Delta_k) + 1
# k = np.linspace(0,kmax,1000000)
k = np.linspace(0,kmax,nkbin)
k[0] = 1E-20
pk_camb = powerspectrum.P_0("~/Work/Thesis/Git/SaclayMocks/etc/PlanckDR12.fits").P(k)
# pk_camb = powerspectrum.P_0("~/Work/Thesis/Git/LyaMocks/bin/data/PlanckDR12.fits").P(k)
r1,xi_camb = powerspectrum.xi_from_pk(k, pk_camb)
p1d_camb = powerspectrum.P_1D(k, pk_camb).P1D(k)

kn = np.pi / d_cell
cut = np.sinc(k*d_cell/(2*np.pi))
r2,xi_camb_cut = powerspectrum.xi_from_pk(k, pk_camb*cut**2)
p1d_camb_cut = powerspectrum.P_1D(k, pk_camb*cut**2).P1D(k)
msk = k < kn
p1d_camb_cut2 = powerspectrum.P_1D(k[msk], pk_camb[msk]).P1D(k)

smooth = np.exp(- k**2 * d_cell**2 / 2)
r3,xi_camb_smooth = powerspectrum.xi_from_pk(k, pk_camb*smooth**2)
r4,xi_camb_cut_smooth = powerspectrum.xi_from_pk(k, pk_camb*(cut**2)*(smooth**2))
p1d_camb_smooth = powerspectrum.P_1D(k, pk_camb*smooth**2).P1D(k)
p1d_camb_cut_smooth = powerspectrum.P_1D(k, pk_camb*cut**2*smooth**2).P1D(k)



# plt.plot(k, pk_camb)
# plt.plot(k, pk_camb*smooth)
# plt.xlabel(r'$k$ [h/Mpc]')
# plt.ylabel('P(k)')
# plt.grid()
# plt.yscale('log')
# plt.xlim(0, kn+0.2)
# plt.ylim(1e-3, 1e5)
# # plt.vlines(kn, ymin=0, ymax=2000, linestyles='--')
# plt.show()

# plt.plot(r1, xi_camb*r1**p, label='camb')
# # plt.plot(r2, xi_camb_cut*r2**p, label='camb+cut')
# plt.plot(r3, xi_camb_smooth*r3**p, label='camb+smoothing')
# # plt.plot(r4, xi_camb_cut_smooth*r4**p, label='camb+cut+smoothing')
# plt.xlabel(r'$r$ [Mpc/h]')
# plt.ylabel(r'$r^2 \xi(r)$')
# plt.grid()
# # plt.legend()
# plt.xlim(-15, 350)
# plt.ylim(-15, 45)
# plt.show()

# f0, ax0 = plt.subplots()
# ax0.plot(k, pk_camb, color='black')
# ax0.plot(k, pk_camb*smooth**2, color='blue')
# ax0.plot(k, pk_camb*cut**2, color='green')
# ax0.plot(k, pk_camb*smooth**2*cut**2, color='red')
# ax0.set_xlabel(r'$k \; [h \; \mathrm{Mpc}^{-1}]$')
# ax0.set_ylabel(r'$P(k) \; [(h^{-1} \mathrm{Mpc})^3]$')
# ax0.grid()
# ax0.set_yscale('log')
# ax0.set_xlim(0, 1.7)
# ax0.set_ylim(1e-3, 1e5)
# plt.tight_layout()

f1, ax1 = plt.subplots()
ax1.plot(k, p1d_camb)
ax1.plot(k, p1d_camb_smooth)
# ax1.plot(k, p1d_camb_cut, color='green')
# ax1.plot(k, p1d_camb_cut2, color='magenta')
# ax1.plot(k, p1d_camb_cut_smooth, color='red')
ax1.set_xlabel(r'$k \; [h \; \mathrm{Mpc}^{-1}]$')
ax1.set_ylabel(r'$P^{\mathrm{1D}}(k) \; [(h^{-1} \mathrm{Mpc})]$')
ax1.grid()
plt.tight_layout()
plt.show()

# f, ax = plt.subplots(nrows=2, sharex=True)
# ax[0].plot(r1, xi_camb*r1**p, label='camb', color='black')
# ax[0].plot(r3, xi_camb_smooth*r3**p, label='camb+smoothing', color='blue')
# # ax[0].plot(r2, xi_camb_cut*r2**p, label='camb+cut', color='green')
# # ax[0].plot(r4, xi_camb_cut_smooth*r4**p, label='camb+cut+smoothing', color='red')

# ax[1].plot(r1, (xi_camb_smooth - xi_camb)*r1**p, color='blue')
# # ax[1].plot(r1, (xi_camb_cut - xi_camb)*r1**p, color='green')
# # ax[1].plot(r1, (xi_camb_cut_smooth - xi_camb)*r1**p, color='red')

# ax[1].set_xlabel(r'$r \;[h^{-1} \mathrm{Mpc}]$')
# ax[0].set_ylabel(r'$r^2 \xi(r)$')
# ax[1].set_ylabel(r'$r^2 (\xi(r) - \xi_{\mathrm{camb}}(r))$')
# # ax[1].set_xlabel(r'$r \;[h^{-1} \mathrm{Mpc}]$')
# # ax[0].set_ylabel(r'$r^2 \xi(r) \; [(h^{-1} \mathrm{Mpc})^2]$')
# # ax[1].set_ylabel(r'$r^2 (\xi(r) - \xi_{\mathrm{camb}}(r)) \; [(h^{-1} \mathrm{Mpc})^2]$')
# ax[0].grid()
# ax[1].grid()
# # ax[0].legend()
# # ax[1].legend()
# plt.xlim(-5,200)
# plt.tight_layout()
plt.show()


