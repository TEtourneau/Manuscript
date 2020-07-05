import numpy as np
import matplotlib.pyplot as plt


z = np.arange(2, 4.3, 0.2)
tau = np.array([0.127, 0.164, 0.203, 0.251, 0.325, 0.386, 0.415, 0.570, 0.716, 0.832, 0.934, 1.061])
tau_err = np.array([0.018, 0.013, 0.009, 0.010, 0.012, 0.014, 0.017, 0.019, 0.030, 0.025, 0.037, 0.091])
trans = np.exp(-tau)
trans_err = trans*tau_err

zz = np.linspace(1.8, 3.6, 1000)
trans_fit = np.exp(-0.00211*(1+zz)**3.7)


plt.errorbar(z, trans, yerr=trans_err, fmt='o', color='black')
plt.plot(zz, trans_fit, linestyle='--', color='blue')
plt.xlabel(r'$z$')
plt.ylabel(r'$\overline{F}$')
plt.grid()
plt.tight_layout()
plt.show()
