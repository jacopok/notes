#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('text.latex', preamble=r'''\usepackage{amsmath}
          \usepackage{physics}
          \usepackage{siunitx}
          ''')
import scipy.constants as sc
from astropy.constants import codata2018 as ac
from astropy.constants import iau2015 as aa
import astropy.units as u
from astropy.cosmology import Planck15 as cosmo
import astropy.uncertainty as aun
          
x_n = np.logspace(-2.5, 2.5)
Q = (ac.m_n - ac.m_p) 
m_n = ac.m_n
m_e = ac.m_e

R_np = x_n ** 3 * (4 * (1 + x_n ** 2) / (x_n ** 4 + 4 * Q / m_n + 4 * (Q ** 2 - m_e ** 2) / m_n ** 2))**(3 / 2)

# plt.loglog(x_n, x_n**3 * (1/R_np  + 1))
# plt.loglog(x_n**3 * (1/R_np  + 1), R_np)
plt.grid()
plt.loglog(x_n, R_np)
plt.xlabel('$x_n$')
plt.ylabel('$R_{\\text{np}}$')
# plt.savefig('neutron-proton-ratio.pdf')
plt.show()
# %%

density_unit = u.g/ u.cm**3

const = (ac.m_p * 8 * np.pi / 3 * (ac.m_p * ac.c / ac.hbar)** 3).to(density_unit)

rho = const * x_n**3 * (1/R_np  + 1)
plt.loglog(rho, R_np)
plt.grid()
plt.xlabel('$\\rho$ (logarithmic, units of g/cm$^3$)')
plt.ylabel('$R_{\\text{np}}$')
plt.savefig('neutron-proton-ratio-density.pdf') 
plt.show()

# %%
