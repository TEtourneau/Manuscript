import numpy as np
from SaclayMocks import powerspectrum
import matplotlib.pyplot as plt 


p = 2
d_cell = 2.19
kmax = 100
k = np.linspace(0,kmax,1000000)
pk_camb = powerspectrum.P_0("~/Work/Thesis/Git/SaclayMocks/etc/PlanckDR12.fits").P(k)
# pk_camb = powerspectrum.P_0("~/Work/Thesis/Git/LyaMocks/bin/data/PlanckDR12.fits").P(k)
r1,xi_camb = powerspectrum.xi_from_pk(k, pk_camb)

kn = np.pi / d_cell
cut = np.sinc(k*d_cell/2)
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

plt.plot(r1, xi_camb*r1**p, label='camb')
plt.plot(r2, xi_camb_cut*r2**p, label='camb+cut')
plt.plot(r3, xi_camb_smooth*r3**p, label='camb+smoothing')
plt.plot(r4, xi_camb_cut_smooth*r4**p, label='camb+cut+smoothing')
plt.xlabel(r'$r$ [Mpc/h]')
plt.ylabel(r'$r^2 \xi(r)$')
plt.grid()
plt.legend()
plt.xlim(-15, 350)
plt.ylim(-15, 45)
plt.show()
