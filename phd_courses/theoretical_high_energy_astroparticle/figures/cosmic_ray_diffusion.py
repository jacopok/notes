import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import codata2018 as ac
import astropy.units as u
from tqdm import tqdm
from scipy.integrate import solve_ivp
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable
from astropy.visualization import quantity_support

gamma = 10
v0 = ac.c * np.sqrt(1 - 1/gamma**2)

# jesus christ who invented Gaussian units
B0 = 2e-6 * u.cm**(-1/2) * u.g**(1/2) / u.s

q = 1 * ac.e.gauss
m = ac.m_p

omega_larmor = (q * B0 / (m * gamma * ac.c)).si

delta_B = B0 / 10000

# set to 1 or -1
sign = 1
# change if needed
phi0 = 0

theta0 = 1.

def derivative_func(k):
    o = omega_larmor.to(1/u.s).value
    kv0 = (sign * k * v0).to(1/u.s).value
    dB_over_B = (delta_B / B0).to(u.dimensionless_unscaled).value
    v_over_c = (v0 / ac.c).to(u.dimensionless_unscaled).value
    def theta_derivative(t, theta):
        return (
            - o
            * dB_over_B
            * v_over_c
            * np.cos(
                    (o - kv0 * np.cos(theta) ) * t 
                    + phi0
                )
            )
    return theta_derivative

    
k_resonance = (sign * omega_larmor / v0 / np.cos(theta0)).to(1/ u.AU)

N_plotted = 20
k_range = np.linspace(.8, 1.2, num=N_plotted) * k_resonance

n_periods = 50
larmor_period = 2 * np.pi / omega_larmor
max_step = (larmor_period / 40).si.value
t_span = (0, (n_periods * larmor_period).si.value)
# t_eval = np.linspace(*t_span, num=50 * n_periods)


def diffusion_over_time():

    cmap = plt.get_cmap('viridis')
    norm = Normalize(min(k_range).value, max(k_range.value))
    mappable = ScalarMappable(norm=norm, cmap=cmap)

        
    for k in tqdm(k_range):
        func = derivative_func(k)
        sol = solve_ivp(func, t_span, y0=[theta0], max_step=max_step)
        
        plt.plot(
            sol.t/larmor_period.si.value, 
            sol.y[0], 
            c=cmap(norm(k.value)), 
            alpha=2**(-1.2*np.log10(N_plotted))
        )
    
    plt.colorbar(mappable=mappable, label=f'k [{k_resonance.unit}]')
    plt.xlabel('Time [Larmor periods]')
    plt.ylabel('Angle [radians]')
    plt.title('Diffusion varying $k$')
    
def final_point_variance():

    theta_final = []
    
    for k in tqdm(k_range):
        func = derivative_func(k)
        sol = solve_ivp(func, t_span, y0=[theta0], 
            max_step=max_step)
        theta_final.append(sol.y[0][-1])

    theta_diffs = (np.array(theta_final) - theta0) * u.rad

    with quantity_support():
        plt.semilogx(k_range, theta_diffs)
        plt.axvline(k_resonance, ls=':', label='Resonance wavenumber')
        plt.grid('on')
    
    plt.legend()
    plt.title(f'$\\Delta \\theta$ over {n_periods} Larmor periods')
    
def omegas():
    with quantity_support():
        plt.semilogx(k_range, omega_larmor - k_range * v0 * np.cos(theta0))
        plt.axvline(k_resonance)
        plt.grid('on')


if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(final_point_variance)
    plot_and_save(diffusion_over_time)
    # plot_and_save(omegas)

    