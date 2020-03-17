import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('text.latex', preamble=r'''\usepackage{amsmath}
          \usepackage{physics}
          \usepackage{siunitx}
          ''')

Omega0 = np.linspace(1.00001, 2)
Omega0_less = np.linspace(.01, .99999)

t = 1 / 2 * Omega0 / (Omega0 - 1)**(3 / 2)\
  * (np.arccos(2 / Omega0 - 1) - 2 / Omega0 * np.sqrt(Omega0 - 1))
t_less = - 1 / 2 * Omega0_less / (1- Omega0_less)**(3 / 2)\
  * (np.arccosh(2 / Omega0_less - 1) - 2 / Omega0_less * np.sqrt(1 - Omega0_less ))

plt.plot(Omega0, t, label='Curved model, $k=1$')
plt.plot(Omega0_less, t_less, label='Curved model, $k=-1$')
plt.axhline(2/3, ls='dashed', lw=1, label='Flat model ')
plt.xlabel('Current ratio of density to critical density $\\Omega_0$')
plt.ylabel('Universe age $t_0$ in units of $H_0^{-1}$')
plt.legend()
plt.savefig('curved_universe_age.pdf', format = 'pdf')
# plt.show()