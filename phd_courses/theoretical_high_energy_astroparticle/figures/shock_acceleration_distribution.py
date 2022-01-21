import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
import matplotlib as mpl

cmap = plt.get_cmap('plasma')

@np.vectorize
def distribution(x, x0, f0, lam):
    if x < 0:
        return max(f0 / (1 - np.exp(x0 / lam)) * (
            np.exp(x / lam)- 
            np.exp(x0 / lam)
        ), 0.)
    else:
        return f0

def shock_acceleration_distribution():
    xs = np.linspace(-5, 5, num=200)
    
    x0 = -3 
    
    f0 = 1. 
    
    # lambda = D / u_1
    lams = [.2, 3.,  10.]
    
    for lam in lams:
        plt.plot(
            xs, 
            distribution(xs, x0, f0, lam), 
            label=f'$D / u_1 = $ {lam}'
            )
    plt.axvline(x=x0, label=f'$x_0 = {x0}$', ls=':', c='black', lw=.5)
    plt.xlabel('$x$')
    plt.ylabel('$f(x)$')
    plt.legend()
    
    
def solve_beta_bvp(k=1., alpha_CR=1., r=4., x0=-1., logp=1.):
    
    spectrum = 3 * r / (r-1)

    def fun(x, y):
        # dy/dx = fun(x, y)
        # y = [beta, beta']
        # beta'' = (beta')**2 - k * beta'
        
        return np.vstack((y[1], y[1] * (y[1] - k) / logp))
    
    def bc(ya, yb):
        # boundary condition residuals
        return np.array([
            ya[0] - alpha_CR, 
            1 - yb[1] * logp / k - yb[0] / spectrum
        ])
    
    xs = np.linspace(x0, 0, num=200)
    beta_guess = k * xs + alpha_CR - k * x0
    beta_guess_derivative = np.gradient(beta_guess, xs)
    ys_guess = np.vstack((beta_guess, beta_guess_derivative))
    
    sol = solve_bvp(fun, bc, x = xs, y = ys_guess, verbose=0, tol=1e-4, max_nodes=100_000)
    print(f'Status = {sol.status}, alpha = {alpha_CR}, logp= {logp}')
    return sol

def cosmic_ray_reacceleration(k, logp):
    r=4.
    x0=-1.
    min_alpha=3
    max_alpha=5
    da = max_alpha - min_alpha
    
    for alpha_CR in np.linspace(min_alpha, max_alpha, num=50):
        
        sol = solve_beta_bvp(r=r, alpha_CR=alpha_CR, x0=x0, k=k, logp=logp)
        plt.plot(sol.x * k, sol.y[0], c=cmap((alpha_CR-min_alpha)/da))
        plt.axhline(3 * r / (r-1), ls=':', c='black', lw=.8)
    plt.xlabel('Rescaled position $X = x u_1 / D$')
    plt.ylabel(r'spectral index $\beta$')
    plt.title(r'$D/u_1= $ ' + str(1/(x0 * k)) + r'$\abs{x_0}$, $\log (p / p_{\text{min}})= $ ' + str(logp))

def beta0_grid(alpha_CR):
    num = 50
    k_range = np.logspace(-1.2, .7, num=num)
    logp_range = np.linspace(.2, 3, num=num)
    
    K, LOGP = np.meshgrid(k_range, logp_range)
    
    beta0 = np.zeros_like(K)
    
    for ind, k in np.ndenumerate(K):
        logp = LOGP[ind]
        print(f'{k=}, {logp=}')
        
        sol = solve_beta_bvp(k=k, logp=logp, alpha_CR=alpha_CR)
        if sol.status == 0:
            beta0[ind] = sol.y[0][-1]
    
    return K, LOGP, beta0    

def show_alpha_spectrum(alpha_CR=3.):
            
    c = plt.contourf(*beta0_grid(alpha_CR), levels=100)
    # plt.clim(alpha_CR, 4)
    plt.xlabel('$k |x_0| = |x_0| u_1 / D$')
    plt.xscale('log')
    plt.ylabel(r'$\log (p / p _{\text{min}})$')
    cbar = plt.colorbar(c, label=r'$\beta (x=0)$: $\alpha_{CR}$ = ' + f'{alpha_CR}')

def show_two_alpha_spectrums():
        
    K, LOGP, beta0_3 = beta0_grid(3.)
    K, LOGP, beta0_5 = beta0_grid(5.)
    
    norm = mpl.colors.Normalize(vmin=3, vmax=5)
    cmap = plt.get_cmap('seismic')
    
    lab = 'varying $p$'
    for k0 in [.5, 1, 2]:
        p = np.logspace(.1, 1)
        k = k0 / np.sqrt(p)
        plt.plot(k, np.log(p), c='black', ls=':', label=lab)
        lab=None
    
    levels = np.linspace(3, 5, num=30)
    
    plt.contour(K, LOGP, beta0_3, cmap=cmap, norm=norm, levels=levels)
    plt.contour(K, LOGP, beta0_5, cmap=cmap, norm=norm, levels=levels)

    plt.legend()
    plt.xlabel('$k |x_0| = |x_0| u_1 / D$')
    plt.xscale('log')
    plt.ylabel(r'$\log (p / p _{\text{min}})$')
    plt.title(r'Final value of $\beta$')
    cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm), label=r'$\beta (x=0)$')


def show_low_alpha_spectrum():
    show_alpha_spectrum(alpha_CR=3.)
def show_high_alpha_spectrum():
    show_alpha_spectrum(alpha_CR=5.)

def cosmic_ray_reacceleration_far():
    cosmic_ray_reacceleration(k=.5, logp=2.)

def cosmic_ray_reacceleration_near():
    cosmic_ray_reacceleration(k=4., logp=.5)


if __name__ == "__main__":
    from make_all_figures import plot_and_save
    # plot_and_save(shock_acceleration_distribution)
    # plot_and_save(cosmic_ray_reacceleration_far)
    # plot_and_save(cosmic_ray_reacceleration_near)
    # plot_and_save(show_low_alpha_spectrum)
    # plot_and_save(show_high_alpha_spectrum)
    plot_and_save(show_two_alpha_spectrums)