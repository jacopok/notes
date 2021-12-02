import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import codata2018 as ac
import astropy.units as u
from tqdm import tqdm
from scipy.integrate import solve_ivp
from matplotlib.colors import Normalize, LogNorm, FuncNorm
from matplotlib.cm import ScalarMappable
from astropy.visualization import quantity_support
from matplotlib.ticker import MultipleLocator
from scipy.integrate import trapezoid
import cmasher as cmr

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
                    (o - kv0 * np.cos(theta)) * t
                + phi0
            )
        )
    return theta_derivative


k_resonance = (sign * omega_larmor / v0 / np.cos(theta0)).to(1 / u.AU)

N_plotted = 200
N_plotted_margins = 100

lim_k_dex = 1

def cubic(x, center=0, width=2*lim_k_dex):
    return (width/2)**(-1.5) * (abs(x - center))**1.5 * np.sign(x-center)
    
def inverse_cubic(y, center=0, width=2*lim_k_dex):
    return center + (abs(y))**(1/1.5) * (width/2) * np.sign(y)

logs = cubic(np.linspace(-lim_k_dex, lim_k_dex, num=N_plotted))

k_range = 10**logs * k_resonance

# k_range = np.concatenate(
#     (
#         np.logspace(-.8, -.3, num=N_plotted_margins),
#         np.logspace(-.3, .3, num=N_plotted),
#         np.logspace(+.3, .8, num=N_plotted_margins),
#     )
# ) * k_resonance


# k_range = np.logspace(-.3, .3, num=N_plotted) * k_resonance

n_periods = 100
larmor_period = 2 * np.pi / omega_larmor

global_pulsations = omega_larmor - k_range * v0 * np.cos(theta0)
integration_oom = 2 * np.pi / (abs(global_pulsations) + omega_larmor)

max_steps = (integration_oom / 4).si.value

t_span = (0, (n_periods * larmor_period).si.value)
t_eval = np.linspace(*t_span, num=4 * n_periods)


def solver():
    for k, dt in tqdm(zip(k_range, max_steps), total=len(k_range)):
        func = derivative_func(k)
        sol = solve_ivp(func, t_span, y0=[theta0], max_step=dt, t_eval=t_eval)
        sol.k = k
        yield sol


def diffusion_over_time():

    cmap = cmr.iceburn
    norm = LogNorm(k_range[0].value, k_range[-1].value)
    # norm = FuncNorm((
    #     lambda x: cubic(np.log10(x/k_resonance.value)),
    #     lambda x: 10**inverse_cubic(x) * k_resonance.value), 
    #     vmin=min(k_range).value, 
    #     vmax=max(k_range).value
    # )

    for sol in solver():
        plt.plot(
            sol.t/larmor_period.si.value,
            sol.y[0],
            c=cmap(norm(sol.k.value)),
            alpha=2**(-1.*np.log10(N_plotted)),
            lw=.5
        )

    mappable = ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(mappable=mappable, label=f'k [{k_resonance.unit}]')
    cbar.ax.hlines(k_resonance.value, 0, 2, color='white')
    # plt.gca().set_facecolor('black')

    plt.xlabel('Time [Larmor periods]')
    plt.ylabel('Angle [radians]')
    plt.title('Diffusion varying $k$')


def final_point_variation():

    theta_final = [sol.y[0][-1] for sol in solver()]

    theta_diffs = (np.array(theta_final) - theta0) * u.rad

    # with quantity_support():
    plt.semilogx(k_range.value, theta_diffs.value, lw=.8)
    plt.axvline(k_resonance.value, ls=':',
                label='Resonance wavenumber', c='black')
    plt.xlabel(f'$k$ [{k_range.unit}]')
    plt.ylabel(f'$\\Delta \\theta$ [{theta_diffs.unit}]')

    plt.grid('on')

    plt.legend()
    plt.title(f'$\\Delta \\theta$ over {n_periods} Larmor periods')


# def final_point_integral():
#     theta_final = [sol.y[0][-1] for sol in solver()]
#     theta_diffs = (np.array(theta_final) - theta0)
#     theta_vars = theta_diffs**2
#     integral = trapezoid(y=theta_vars, x=k_range.value)
#     th_value = (omega_larmor**2 * (delta_B / B)**2 * np.pi / v0 / np.cos(theta0) * t_span[1]).si
#     print(th_value)
#     print(integral)
    
    # # with quantity_support():
    # plt.semilogx(k_range.value, theta_vars, lw=.8)
    # plt.axvline(k_resonance.value, ls=':',
    #             label='Resonance wavenumber', c='black')
    # plt.xlabel(f'$k$ [{k_range.unit}]')
    # plt.ylabel(f'$\\Delta \\theta$ [{theta_diffs.unit}]')

    # plt.grid('on')

    # plt.legend()
    # plt.title(f'$\\Delta \\theta$ over {n_periods} Larmor periods')



def integration_periods_plot():
    with quantity_support():
        plt.semilogx(k_range, integration_oom)
        plt.axvline(k_resonance)
        plt.grid('on')


if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(final_point_variation)
    # plot_and_save(diffusion_over_time)
    # plot_and_save(integration_periods_plot)
    # final_point_integral()