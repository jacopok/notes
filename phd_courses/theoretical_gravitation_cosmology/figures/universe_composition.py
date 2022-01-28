import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.constants as sc
from astropy.constants import codata2018 as ac
from astropy.constants import iau2015 as aa
import astropy.units as u
from astropy.cosmology import Planck15 as cosmo
import astropy.uncertainty as aun

def universe_composition():
    z = np.concatenate(([0], np.logspace(-5, 6, num=1000)))
    
    components = {
        'Dark energy': cosmo.Ode(z),
        'Dark matter': cosmo.Odm(z),
        'Baryonic matter': cosmo.Ob(z),
        'Radiation': cosmo.Ogamma(z) + cosmo.Onu(z)
        # 'Photon radiation': cosmo.Ogamma(z),
        # 'Neutrino radiation': cosmo.Onu(z)
    }
    colors = {
        'Dark energy': 'black',
        'Dark matter': '#361be3',
        'Baryonic matter': 'red',
        'Radiation': 'yellow'
        # 'Photon radiation': 'yellow',
        # 'Neutrino radiation': 'orange',
    }
    alphas = {
        'Dark energy': 1.,
        'Dark matter': 1.,
        'Baryonic matter': .5,
        'Radiation': .5,
        # 'Photon radiation': .5,
        # 'Neutrino radiation': .5,
    }
    
    current_budget = np.ones_like(z)
    
    for name, comp in components.items():
        
        plt.fill_between(z,
            y1=current_budget-comp,
            y2=current_budget,
            alpha=alphas[name],
            label=name, 
            color=colors[name])
        current_budget -= comp

    plt.xlabel('Redshift $z$')
    plt.ylabel('Energy budget: decomposition of $\Omega = \sum_i \Omega_i$')
    plt.ylim(0, 1)
    plt.xlim(min(z), max(z))
    plt.legend(loc='upper left')
    
    ax = plt.gca()
    linthresh = 1e-3
    ax.set_xscale('symlog', linthresh=1e-3)
    ax.xaxis.set_minor_locator(matplotlib.ticker.LogLocator(base = 10.0, subs = np.arange(.1, 1., .1), numticks = 12))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    
    new_ticks = [tick for tick in ax.get_xticks(minor=True) if tick >= linthresh/10 and tick <= max(z)]
    ax.set_xticks(ticks = new_ticks, minor=True)

def hubble_rate():
    z = np.logspace(-2, 6)
    
    plt.loglog(z, cosmo.H(z))

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(universe_composition)
    plot_and_save(hubble_rate)