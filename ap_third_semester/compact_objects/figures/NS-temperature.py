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

theta = np.linspace(0, np.pi / 2)
cos_big_theta = 2 * np.cos(theta) / np.sqrt(1 + 3 * np.cos(theta)** 2)

plt.plot(theta, np.sqrt(np.abs(cos_big_theta)))
plt.xlabel('$\\theta$')
plt.ylabel('$T / T_P$')
plt.savefig('NS-temperature.pdf')