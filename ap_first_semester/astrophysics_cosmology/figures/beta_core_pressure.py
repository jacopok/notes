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

a = np.pi ** 2 * ac.k_B**2 / 15 / ac.hbar**3 / ac.c**3

beta = np.linspace(0, 1, num=200)
mbar = .62

M = np.sqrt(ac.G ** (-3 ) * (np.pi / 36)**(-1) * (3/a * (1-beta) / beta**4) * (ac.k_B / mbar)**(4))

with quantity_support():
    plt.semilogx(M.decompose(), beta)
    plt.show()
    
    # the units don't work!