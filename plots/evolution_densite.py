import numpy as np
import matplotlib.pyplot as plt

SMALL_SIZE = 18
MEDIUM_SIZE = 20
BIGGER_SIZE = 22
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rc('axes', labelsize=MEDIUM_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)
plt.rc('figure', figsize=(9,7))

zeq = 3411
zstar = 1090
zdrag = 1059

a = np.logspace(-8,0,1000)
omega_m = 0.3111
omega_l = 0.6889
omega_r = omega_m / (1+zeq)
rad = omega_r * a**-4
mat = omega_m * a**-3
cc = omega_l * np.ones_like(a)

# z = np.logspace(-10, 10, 1000000)
# omega_m = 0.3111
# omega_l = 0.6889
# omega_r = 1e-4
# rad = omega_r * (1+z)**4
# mat = omega_m * (1+z)**3
# cc = omega_l * np.ones_like(z)


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

ax1.plot(a, rad)
ax1.plot(a, mat)
ax1.plot(a, cc)
ax1.set_xlabel('a')
# ax1.plot(z, rad)
# ax1.plot(z, mat)
# ax1.plot(z, cc)
# ax1.xlim(1e10, 1)
# ax1.xlabel('z')

ax1.grid()
ax1.set_xscale('log')
ax1.set_yscale('log')

new_tick_locations = np.array([1e-5, 1e-4, 1/(1+zeq), 1/(1+zstar), 1/(1+0.3), 1])
new_names = ['1e5', '1e4', r'z_{eq}', r'z_{s}', '0.3', '0']

def tick_function(X):
    V = 1/X - 1
    return ["%.3f" % z for z in V]

ax2.set_xlim(ax1.get_xlim())
ax2.set_xlabel("z")
ax2.set_xscale('log')
ax2.set_xticks(new_tick_locations)
# ax2.set_xticklabels(tick_function(new_tick_locations))
ax2.set_xticklabels(new_names)

plt.tight_layout()
plt.show()

# plt.plot(a, rad)
# plt.plot(a, mat)
# plt.plot(a, cc)
# plt.xlabel('a')
# # plt.plot(z, rad)
# # plt.plot(z, mat)
# # plt.plot(z, cc)
# # plt.xlim(1e10, 1)
# # plt.xlabel('z')

# plt.grid()
# plt.ylabel(r'$\rho \,\mathrm{[a. u.]}$')
# plt.xscale('log')
# plt.yscale('log')
# plt.tight_layout()
# plt.show()
