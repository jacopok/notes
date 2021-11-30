import numpy as np
import nptyping as npt
from typing import Any, Optional
import matplotlib.pyplot as plt
from astropy.constants import codata2018 as ac
import astropy.units as u
from astropy.visualization import quantity_support

neutrino_names = ['e', 'mu', 'tau']

def osc_probability(
        E_over_L: np.float64,
        initial_flavour: int = 0,
        final_flavour: int = 1,
        masses: Optional[npt.NDArray[(3,), np.float64]] = None,
        mixing_matrix: Optional[npt.NDArray[(3, 3), np.complex128]] = None,
    ):
    """
    Probability of oscillation

    Args:
        E_over_L (np.float64): in eV**2
        masses (npt.NDArray): [m1, m2, m3], in eV
        mixing_matrix (npt.NDArray): [flavour, mass], where flavour is in (e, mu, tau)
        initial_flavour (int): number between 0=e, 1=mu, 2=tau
        final_flavour (int): number between 0=e, 1=mu, 2=tau
    """
    
    if masses is None and mixing_matrix is None:
        masses, mixing_matrix = prepare_oscillation_parameters()

    return abs(np.sum(
        np.conj(mixing_matrix[initial_flavour, :, np.newaxis]) *
        mixing_matrix[final_flavour, :, np.newaxis] *
        np.exp(- 1j * masses[:, np.newaxis]**2 
                    / E_over_L[np.newaxis, :] 
                    / 2),
        axis=0
    ))**2


def matrix(theta_12, theta_13, theta_23, delta_CP):

    c12 = np.cos(theta_12)
    c13 = np.cos(theta_13)
    c23 = np.cos(theta_23)

    s12 = np.sin(theta_12)
    s13 = np.sin(theta_13)
    s23 = np.sin(theta_23)

    return np.array([
        [
            c12*c13, 
            s12*c13, 
            s13*np.exp(-1j*delta_CP)
        ],
        [
            -s12*c23-c12*s23*s13*np.exp(1j*delta_CP),
            c12*c23 - s12*s13*s23*np.exp(1j * delta_CP), 
            s23*c13
            ],
        [
            s12*s23 - c12*c23*s13*np.exp(1j*delta_CP), 
            -c12*s23 - s12*c23*s13*np.exp(1j * delta_CP),
            c23*c13
        ]
    ])

def osc_plot(E_over_L, op, initial, final):

    with quantity_support():
        
        plt.semilogx(E_over_L, op, label=f'Transition from {neutrino_names[initial]} to {neutrino_names[final]}')

        # ax = plt.gca()
        # ax2 = ax.twiny()
        # ax2.semilogx(1 / E_over_L, op, alpha=.1)
        # ax2.set_xlim(reversed(ax2.get_xlim()))

def flux_plot(L_over_E, op, flux):

    with quantity_support():
        plt.plot(L_over_E, (1 - op) * flux) 
        plt.plot(L_over_E, flux, ls=':')

def flux(E):
    """ Parametrization from https://arxiv.org/abs/0807.3203, 
    where the energy is in MeV (?)
    """
    
    # energy and momentum of the emitted positron
    Ep = E - 1.29333236
    pp = np.sqrt(Ep**2 - 0.51099895**2)
    
    return (
        .58 * np.exp(.87 - .16*E - .091*E**2)+
        .30 * np.exp(.896 - .239*E - .0981*E**2)+
        .07 * np.exp(.976 - .162*E - .079*E**2)+
        .05 * np.exp(.793 - .08*E - .1085*E**2)
        ) * Ep * pp

def prepare_oscillation_parameters():

    # in terms of the given values for sin**2 (theta)
    theta_12 = np.arcsin(np.sqrt(.31))
    theta_13 = np.arcsin(np.sqrt(.558))
    theta_23 = np.arcsin(np.sqrt(.02241))
    delta_CP = (222 * u.degree).to(u.rad).value
    
    U = matrix(theta_12, theta_13, theta_23, delta_CP)
    
    m1 =  (2 * u.meV).to(u.eV)
    
    delta_m21_square = (73.9 * u.meV**2).to(u.eV**2)
    delta_m32_square = (2449 * u.meV**2).to(u.eV**2)
    
    m2 = np.sqrt(delta_m21_square - m1**2)
    m3 = np.sqrt(delta_m32_square - m2**2)

    masses = np.array([m1.value, m2.value, m3.value])

    return masses, U

def oscillations():

    # parameters for Juno
    E_range = np.logspace(0, 4, num=2000) * u.MeV
    L = 50 * u.km
    
    E_over_L_range = (E_range / L * ac.c * ac.hbar).to(u.eV**2).value

    initial = 0
    rescaled_E_over_L = (E_over_L_range * u.eV**2 / ac.hbar / ac.c).to(u.MeV / u.km)

    for final in range(3):
        op = osc_probability(E_over_L_range, initial, final)
        osc_plot(rescaled_E_over_L, op, initial, final)
    
    plt.legend()

def juno_flux():

    E = 1 * u.MeV 
    L_range = np.linspace(5, 32, num=1000) * u.km
    
    E_over_L_range = (E / L_range * ac.c * ac.hbar).to(u.eV**2).value
    
    op = osc_probability(E_over_L_range, 0, 0)
    
    rescaled_E_over_L = (E_over_L_range * u.eV**2 / ac.hbar / ac.c).to(u.MeV / u.km)
    
    f = flux((rescaled_E_over_L * 60 * u.km).to(u.MeV).value)
    
    op = np.convolve(op, np.ones(25)/25, mode='same')

    flux_plot(1/rescaled_E_over_L, op, f)
    plt.ylabel('Flux [arbitrary units]')


if __name__ == "__main__":

    from make_all_figures import plot_and_save
    plot_and_save(oscillations)
    plot_and_save(juno_flux)
    