import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('text.latex', preamble=r'''\usepackage{amsmath}
          \usepackage{physics}
          \usepackage{siunitx}
          ''')
from matplotlib import ticker
from tqdm import tqdm

import scipy.constants as sc
from astropy.constants import codata2018 as ac
from astropy.constants import iau2015 as aa
import astropy.units as u
from astropy.cosmology import Planck15 as cosmo

from astropy.visualization import quantity_support

x = np.linspace(0,10, num=200)

I = 3 / 2 / x ** 4 * (x * (1 + x ** 2)**(1 / 2)* (2 * x ** 2 / 3 - 1) + np.log(x + (1 + x ** 2)**(1 / 2)))

plt.plot(I ** (3 / 2), x)
plt.ylabel('$x_F$')
plt.xlabel('$M / M _{\\text{Ch}} = I(x_F)^{3/2}$')
plt.savefig('chandrasekhar_limit.pdf')
