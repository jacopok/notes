import numpy as np  # type: ignore
import nptyping as npt # type: ignore
from typing import Any, Optional, Tuple
import matplotlib.pyplot as plt # type: ignore
from astropy.constants import codata2018 as ac # type: ignore
import astropy.units as u # type: ignore
from astropy.visualization import quantity_support  # type: ignore

neutrino_names = ['e', 'mu', 'tau']

m_min = 0 * u.eV

def osc_probability(
    E_over_L: npt.NDArray[float],
    initial_flavour: int = 0,
    final_flavour: int = 1,
    ordering: str = 'normal',
    ):
    """
    Probability of oscillation

    Args:
        E_over_L (float): in eV**2
        initial_flavour (int): number between 0=e, 1=mu, 2=tau
        final_flavour (int): number between 0=e, 1=mu, 2=tau
        ordering (str): either normal or inverted
    """

    masses, mixing_matrix = prepare_oscillation_parameters(ordering)

    return abs(np.sum(
        np.conj(mixing_matrix[initial_flavour, :, np.newaxis]) *
        mixing_matrix[final_flavour, :, np.newaxis] *
        np.exp(- 1j * masses[:, np.newaxis]**2
               / E_over_L[np.newaxis, :]
               / 2),
        axis=0
    ))**2


def matrix(theta_12: float, theta_13: float, theta_23: float, delta_CP: float):

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


@u.quantity_input
def osc_plot(
    E_over_L,
    # E_over_L: u.Quantity[u.MeV / u.km],
    op: np.ndarray, 
    initial: int, 
    final: int
    ):

    with quantity_support():

        plt.semilogx(
            E_over_L, op, label=f'Transition from {neutrino_names[initial]} to {neutrino_names[final]}')

@u.quantity_input
def flux_plot(
    L_over_E,
    # L_over_E: u.Quantity[u.km / u.MeV],
    op: np.ndarray, 
    flux: np.ndarray, 
    **kwargs):

    with quantity_support():
        plt.plot(L_over_E, op * flux, **kwargs)
        plt.plot(L_over_E, flux, ls=':')


def flux(E: np.ndarray):
    """ Parametrization from https://arxiv.org/abs/0807.3203, 
    where the energy is in MeV (?)
    """
    # E /= 1.2

    # energy and momentum of the emitted positron
    Ep = E - 1.29333236  # Energy minus neutron-proton mass difference
    pp = np.sqrt(Ep**2 - 0.51099895**2)  # proton energy minus positron mass

    return (
        .58 * np.exp(.87 - .16*E - .091*E**2) +
        .30 * np.exp(.896 - .239*E - .0981*E**2) +
        .07 * np.exp(.976 - .162*E - .079*E**2) +
        .05 * np.exp(.793 - .08*E - .1085*E**2)
    ) * Ep * pp * E**2


def prepare_oscillation_parameters(ordering: str = 'normal') -> Tuple[npt.NDArray[(3, ), float], npt.NDArray[(3, 3), np.complex128]]:

    # based on the values given in http://arxiv.org/abs/1507.05613

    if ordering == 'normal':
        theta_12 = np.arcsin(np.sqrt(.308))
        theta_13 =  np.arcsin(np.sqrt(.1)) / 2
        # theta_13 = np.arcsin(np.sqrt(.0234))
        theta_23 = np.arcsin(np.sqrt(.437))
        delta_m21_square = 7.54e-5 * u.eV**2
        delta_m31_square = 2.47e-3 * u.eV**2
        delta_CP = 1.39 * np.pi

        m1 = m_min
        m2 = np.sqrt(delta_m21_square + m1**2)
        m3 = np.sqrt(delta_m31_square + m1**2)

    elif ordering == 'inverted':
        theta_12 = np.arcsin(np.sqrt(.308))
        theta_13 = np.arcsin(np.sqrt(.0240))
        theta_23 = np.arcsin(np.sqrt(.455))
        delta_m21_square = 7.54e-5 * u.eV**2
        delta_m31_square = -2.42e-3 * u.eV**2
        delta_CP = 1.31 * np.pi

        m3 = m_min
        m1 = np.sqrt(m3**2 - delta_m31_square)
        m2 = np.sqrt(delta_m21_square + m1**2)

    else:
        raise NotImplementedError(f'Ordering {ordering} not found')

    U = matrix(theta_12, theta_13, theta_23, delta_CP)
    masses = np.array([m1.value, m2.value, m3.value])

    return masses, U


def oscillations():

    # parameters for Juno
    E_range = np.logspace(0, 4, num=2000) * u.MeV
    L = 60 * u.km

    E_over_L_range = (E_range / L * ac.c * ac.hbar).to(u.eV**2).value

    initial = 0
    rescaled_E_over_L = (E_over_L_range * u.eV**2 /
                         ac.hbar / ac.c).to(u.MeV / u.km)

    for final in range(3):
        op = osc_probability(E_over_L_range, initial, final)
        osc_plot(rescaled_E_over_L, op, initial, final)

    plt.legend()


def juno_flux():

    E = 1 * u.MeV
    L_range = np.linspace(5, 32, num=1000) * u.km

    E_over_L_range = (E / L_range * ac.c * ac.hbar).to(u.eV**2).value

    op_normal = osc_probability(E_over_L_range, 0, 0, ordering='normal')
    op_inverted = osc_probability(E_over_L_range, 0, 0, ordering='inverted')

    rescaled_E_over_L = (E_over_L_range * u.eV**2 /
                         ac.hbar / ac.c).to(u.MeV / u.km)

    f = flux((rescaled_E_over_L * 60 * u.km).to(u.MeV).value)

    # op = np.convolve(op, np.ones(25)/25, mode='same')

    flux_plot(1/rescaled_E_over_L, op_normal, f, label='normal')
    flux_plot(1/rescaled_E_over_L, op_inverted, f, label='inverted')
    plt.ylabel('Flux [arbitrary units]')
    plt.legend()

if __name__ == "__main__":

    from make_all_figures import plot_and_save
    plot_and_save(oscillations)
    plot_and_save(juno_flux)
