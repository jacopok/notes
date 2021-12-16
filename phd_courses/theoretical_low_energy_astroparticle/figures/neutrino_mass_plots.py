import numpy as np
import matplotlib.pyplot as plt

import scipy.constants as sc
from astropy.constants import codata2018 as ac
from astropy.constants import iau2015 as aa
import astropy.units as u
from astropy.cosmology import Planck15 as cosmo
import astropy.uncertainty as aun
from astropy.visualization import quantity_support

m0 = np.logspace(-4, 0) * u.eV

small_dm = 7.5e-5 * u.eV**2
big_dm = 2.5e-3 * u.eV**2

def mass_plot_no():
    
    m2 = np.sqrt(m0**2 + small_dm)
    m3 = np.sqrt(m0**2 + big_dm)
    
    with quantity_support():
        
        plt.loglog(m0, m0)
        plt.loglog(m0, m2)
        plt.loglog(m0, m3)
        plt.loglog(m0, np.sqrt(.7 * m0**2 + .3 * m2**2 + .02 * m3**2))
    
def mass_plot_io():
    
    m1 = np.sqrt(m0**2 + big_dm)
    m2 = np.sqrt(m1**2 + small_dm)
    
    with quantity_support():

        plt.loglog(m0, m0)
        plt.loglog(m0, m1)
        plt.loglog(m0, m2)
        plt.loglog(m0, np.sqrt(.7 * m1**2 + .3 * m2**2 + .02 * m0**2))

def majorana_mass():
    
    pass

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(mass_plot_no)
    plot_and_save(mass_plot_io)