#%%

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

log_T = np.linspace(4, 15)

urca = 27 + 6 * (log_T - 9)
murca = 21 + 8 * (log_T - 9)

plt.plot(log_T, urca, label='URCA process')
plt.plot(log_T, murca, label='MURCA process')
plt.xlabel('$\\log_{10} T$ [K]')
plt.ylabel('$\\log_{10} \\epsilon_\\nu $ [erg / cm$^3$ / s]')
plt.legend()
plt.savefig('urca-murca.pdf')