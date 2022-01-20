import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

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
    xs = np.linspace(-5, 5, num=1000)
    
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
    
    
def solve_beta_bvp(k=1., alpha_CR=1., r=4., x0=-1.):
    
    spectrum = 3 * r / (r-1)

    def fun(x, y):
        # dy/dx = fun(x, y)
        # y = [beta, beta']
        # beta'' = (beta')**2 - k * beta'
        
        return np.vstack((y[1], y[1] * (y[1] - k)))
    
    def bc(ya, yb):
        # boundary condition residuals
        return np.array([
            ya[0] - alpha_CR, 
            1 - yb[1] / k - yb[0] / spectrum
        ])
    
    xs = np.linspace(x0, 0, num=200)
    beta_guess = k * xs + alpha_CR - k * x0
    beta_guess_derivative = np.gradient(beta_guess, xs)
    ys_guess = np.vstack((beta_guess, beta_guess_derivative))
    
    sol = solve_bvp(fun, bc, x = xs, y = ys_guess, verbose=2, tol=1e-4)
    return sol

def cosmic_ray_reacceleration(k):
    r=4.
    x0=-1.
    for alpha_CR in np.linspace(3, 5, num=50):
        
        sol = solve_beta_bvp(r=r, alpha_CR=alpha_CR, x0=x0, k=k)
        plt.plot(sol.x * k, sol.y[0], c=cmap((alpha_CR-3)/2))
        plt.axhline(3 * r / (r-1), ls=':', c='black', lw=.8)
    plt.xlabel('Rescaled position $X = x u_1 / D$')
    plt.ylabel(r'spectral index $\beta$')
    plt.title(f'$X_0$ = {x0 * k}')

def cosmic_ray_reacceleration_far():
    cosmic_ray_reacceleration(5.)

def cosmic_ray_reacceleration_near():
    cosmic_ray_reacceleration(.2)


if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(shock_acceleration_distribution)
    plot_and_save(cosmic_ray_reacceleration_far)
    plot_and_save(cosmic_ray_reacceleration_near)
    