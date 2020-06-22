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

kn = np.pi / d_cell
cut = np.sinc(k*d_cell/(2*np.pi))
r2,xi_camb_cut = powerspectrum.xi_from_pk(k, pk_camb*cut**2)

smooth = np.exp(- k**2 * d_cell**2 / 2)
r3,xi_camb_smooth = powerspectrum.xi_from_pk(k, pk_camb*smooth**2)
r4,xi_camb_cut_smooth = powerspectrum.xi_from_pk(k, pk_camb*(cut**2)*(smooth**2))

# plt.plot(k, pk_camb)
# plt.plot(k, pk_camb*smooth)
# plt.xlabel(r'$k$ [h/Mpc]')
# plt.ylabel('P(k)')
# plt.grid()
# plt.yscale('log')
# plt.xlim(-0.2, kn+0.2)
# plt.vlines(kn, ymin=0, ymax=2000, linestyles='--')
# plt.show()

# plt.plot(r1, xi_camb*r1**p, label='camb')
# plt.plot(r2, xi_camb_cut*r2**p, label='camb+cut')
# plt.plot(r3, xi_camb_smooth*r3**p, label='camb+smoothing')
# plt.plot(r4, xi_camb_cut_smooth*r4**p, label='camb+cut+smoothing')
# plt.xlabel(r'$r$ [Mpc/h]')
# plt.ylabel(r'$r^2 \xi(r)$')
# plt.grid()
# plt.legend()
# plt.xlim(-15, 350)
# plt.ylim(-15, 45)
# plt.show()

f, ax = plt.subplots(nrows=2, sharex=True)
ax[0].plot(r1, xi_camb*r1**p, label='camb', color='black')
ax[0].plot(r3, xi_camb_smooth*r3**p, label='camb+smoothing', color='blue')
ax[0].plot(r4, xi_camb_cut_smooth*r4**p, label='camb+cut+smoothing', color='red')

ax[1].plot(r1, (xi_camb_smooth - xi_camb)*r1**p, color='blue')
ax[1].plot(r1, (xi_camb_cut_smooth - xi_camb)*r1**p, color='red')
ax[1].plot(r1, (xi_camb_cut_smooth - xi_camb_smooth)*r1**p, color='green')

ax[1].set_xlabel(r'$r$ [Mpc/h]')
ax[0].set_ylabel(r'$r^2 \xi(r)$')
ax[1].set_ylabel(r'$r^2 (\xi(r) - \xi_{\mathrm{camb}}(r))$')
ax[0].grid()
ax[1].grid()
# ax[0].legend()
# ax[1].legend()
plt.xlim(-5,200)
plt.tight_layout()
plt.show()
