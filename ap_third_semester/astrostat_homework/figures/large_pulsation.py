import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('text.latex', preamble=r'''\usepackage{amsmath}
          \usepackage{physics}
          \usepackage{siunitx}
          ''')

N = 500
ts = np.arange(N)
omegas = np.logspace(-5, 5, num=10000)
h = np.average(np.sin(2 * omegas[:, np.newaxis] * ts[np.newaxis,:]), axis=1)
c = np.average(np.cos(omegas[:, np.newaxis] * ts[np.newaxis,:])**2, axis=1)
s = np.average(np.sin(omegas[:, np.newaxis] * ts[np.newaxis,:])**2, axis=1)
plt.semilogx(omegas, c, lw=.7, label='$c/N$')
plt.semilogx(omegas, h, lw=.7, label='$h/N$')
plt.semilogx(omegas, s, lw=.7, label='$s/N$')
plt.xlabel('$\\omega$')
plt.grid('on')
# plt.axvline(1, label='$\\Delta^{-1}$', c='black', ls=':')
plt.axvline(2 * np.pi / 2, label='$\\pi \\Delta^{-1}$ = Nyquist pulsation', c='black', ls=':', lw=.7)
plt.axvline(1/N, c='black', ls='--', label='$1/(N\\Delta )$', lw=.7)
# plt.axvspan(1/np.pi, np.pi, color='blue', alpha=.1, label='Approximation works well')
plt.legend()
plt.savefig('large_pulsation.pdf')
plt.close()
