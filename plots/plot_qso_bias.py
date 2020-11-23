import numpy as np
import matplotlib.pyplot as plt

z = np.array([1.06, 1.35, 1.65, 1.99])
bias = np.array([1.75, 2.06, 2.57, 3.03])
bias_err = np.array([0.08, 0.08, 0.09, 0.11])

z_mean = 1.55
bias_mean = 2.43
bias_err_mean = 0.05

z_param = np.linspace(1,4)
param1 = 3.7 * ((1+z_param)/(1+2.33))**1.7
param2 = 3.58 * ((1+z_param)/(1+2.33))**1.52
param3 = 3.77 * ((1+z_param)/(1+2.334))**1.44

plt.errorbar(z, bias, yerr=bias_err, fmt='.', label='Laurent+2017', color='black')
plt.errorbar(z_mean, bias_mean, yerr=bias_err_mean, fmt='x', label='Laurent+2017', color='black')
plt.plot(z_param, param1, label='parametrisation1')
plt.plot(z_param, param2, label='parametrisation2')
plt.plot(z_param, param3, label='parametrisation3')
plt.grid()
plt.legend()
plt.xlabel('z')
plt.ylabel(r'$b_{\mathrm{QSO}}$')
plt.show()
