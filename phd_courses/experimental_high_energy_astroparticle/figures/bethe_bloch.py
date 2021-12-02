import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import codata2018 as ac
import astropy.units as u

Z = 10
A = 20
rho = 1 * u.g / u.cm**3

n = (ac.N_A * Z * rho / A / ac.u).si
z = 1


def bethe_bloch_func(beta):
    gamma = 1/ np.sqrt(1 - beta**2)
    
    prefactor = (
        4 * np.pi 
        / (ac.m_e * ac.c**2)
        * n * z**2
        / beta**2
        * (ac.e**2
            / (4 * np.pi * ac.eps0)
        )**2
    ).si
    
    return prefactor * 

def bethe_bloch():
    beta = np.logspace(.1, 3)
    

def energy_loss():
    pass

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(bethe_bloch)
    plot_and_save(energy_loss)