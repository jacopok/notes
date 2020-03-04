import numpy as np
import scipy.constants as sc
import astropy.constants as ac
import astropy.units as u
from astropy.cosmology import Planck15 as cosmo
import astropy.uncertainty as aun
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

R = ac.k_B / ac.u
mu = .62
T=  1e7*u.K

a = np.sqrt(R * T / mu)

fs = np.logspace(-5, 0, num=1000) * ac.GM_sun / ac.R_sun**2

def r(f):
  numerator = np.sqrt(1 + f * ac.GM_sun / a ** 4) - 1
  denominator = f / a ** 2
  return (numerator / denominator)
  
plt.plot(fs[1:] / ac.GM_sun * ac.R_sun ** 2, r(fs[1:]) / ac.R_sun)
plt.xlabel("Normalized force: $f / (GM_\\odot / R_\\odot^2)$")
plt.ylabel("Normalized critical radius: $r_c / R_\\odot$")
plt.show(block=False)
