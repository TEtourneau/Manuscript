import matplotlib.pyplot as plt
import numpy as np


t = np.linspace(0, 2.5, 100000)

f, axs = plt.subplots(3,2, sharex=True)

f1 = 1
f2 = 2.2
f3 = 0.3
a1 = 1
a2 = 0.5
a3 = 2

axs[0, 0].plot(t, a1*np.sin(f1*2*np.pi*t), color='r')
axs[0, 1].bar(f1, a1, width=0.05, color='b')

axs[1, 0].plot(t, a1*np.sin(f1*2*np.pi*t), color='mediumblue', linestyle='--')
axs[1, 0].plot(t, a2*np.sin(f2*2*np.pi*t), color='green', linestyle='--')
axs[1, 0].plot(t, a1*np.sin(f1*2*np.pi*t) + a2*np.sin(f2*2*np.pi*t), color='r')
axs[1, 1].bar(f1, a1, width=0.05, color='b')
axs[1, 1].bar(f2, a2, width=0.05, color='b')

axs[2, 0].plot(t, a1*np.sin(f1*2*np.pi*t), color='mediumblue', linestyle='--')
axs[2, 0].plot(t, a2*np.sin(f2*2*np.pi*t), color='green', linestyle='--')
axs[2, 0].plot(t, a3*np.sin(f3*2*np.pi*t), color='c', linestyle='--')
axs[2, 0].plot(t, a1*np.sin(f1*2*np.pi*t) + a2*np.sin(f2*2*np.pi*t) + a3*np.sin(f3*2*np.pi*t), color='r')
axs[2, 1].bar(f1, a1, width=0.05, color='b')
axs[2, 1].bar(f2, a2, width=0.05, color='b')
axs[2, 1].bar(f3, a3, width=0.05, color='b')

axs[0, 1].set_xlim(0, t.max())
axs[1, 1].set_xlim(0, t.max())
axs[2, 1].set_xlim(0, t.max())

axs[0, 0].grid()
axs[1, 0].grid()
axs[2, 0].grid()
axs[0, 0].set_title('Espace temporel')
axs[0, 1].set_title('Espace fréquenciel')
# axs[0, 0].set_xlabel(r'$t [s]$')
# axs[1, 0].set_xlabel(r'$t [s]$')
# axs[0, 1].set_xlabel(r'$f [Hz]$')
# axs[1, 1].set_xlabel(r'$f [Hz]$')
axs[2, 0].set_xlabel(r'$\mathrm{t}\, [\mathrm{s}]$')
axs[2, 1].set_xlabel(r'$\mathrm{f}\, [\mathrm{Hz}]$')

axs[0, 0].set_ylabel(r'$\mathrm{y}(\mathrm{t})$')
axs[1, 0].set_ylabel(r'$\mathrm{y}(\mathrm{t})$')
axs[2, 0].set_ylabel(r'$\mathrm{y}(\mathrm{t})$')
axs[0, 1].set_ylabel(r'$\tilde{\mathrm{y}}(\mathrm{f})$')
axs[1, 1].set_ylabel(r'$\tilde{\mathrm{y}}(\mathrm{f})$')
axs[2, 1].set_ylabel(r'$\tilde{\mathrm{y}}(\mathrm{f})$')

plt.tight_layout()
plt.show()
