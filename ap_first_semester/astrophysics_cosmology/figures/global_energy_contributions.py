import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import scipy.constants as sc
from astropy.constants import codata2018 as ac
from astropy.constants import iau2015 as aa
import astropy.units as u
from astropy.cosmology import Planck15 as cosmo
import astropy.uncertainty as aun

a = np.linspace(.01, 3, num=1000)
# a=1 now

Om0 = cosmo.Om0
Olambda0 = 1 - cosmo.Om0 - cosmo.Ob0
Orad0 = (2.5e-5 * u.littleh**-2).to(u.dimensionless_unscaled, equivalencies=u.with_H0())

def get_rho(a, w, O0):
  return(O0 * a**(-3*(1+w)))

plt.plot(a, get_rho(a, -1, Olambda0), label='Cosmological constant')
plt.plot(a, get_rho(a, 1/3, Orad0), label='Radiation')
plt.plot(a, get_rho(a, 0, Om0), label='Matter')
plt.xlabel('Scale factor: $a / a_0$')
plt.ylabel('Energy density: $\\rho / \\rho_{0c}$')
plt.axvline(x=1, label='Now', c='black')
plt.ylim(0,2)
plt.legend()
plt.savefig('global_energy_contributions.pdf', format = 'pdf')