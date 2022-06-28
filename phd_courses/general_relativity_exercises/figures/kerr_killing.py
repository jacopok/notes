import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from tqdm import tqdm

M = 1
a = .99
cmap = plt.get_cmap('plasma')

def kerr_killing():
    
    thetas = np.linspace(1e-2, np.pi/2, num=10)
    minalpha = 1
    maxalpha = 0
    for theta in tqdm(thetas):
        for r in r_range(theta, num=int((theta + .1) * 50)):
            alpha = best_alpha(r, theta)
            minalpha = min(alpha, minalpha)
            maxalpha = max(alpha, maxalpha)


    thetas = np.linspace(1e-2, np.pi/2, num=250)
    norm = Normalize(minalpha, maxalpha)
    for theta in tqdm(thetas):
        for r in r_range(theta, num=int((theta + .1) * 60)):
            alpha = best_alpha(r, theta)
            plt.scatter(
                r*np.sin(theta), 
                r*np.cos(theta), 
                color=cmap(norm(alpha)),
                s=1.2
            )

    plt.gca().set_aspect('equal')

    plt.colorbar(ScalarMappable(cmap=cmap, norm=norm), label=r'$\alpha$')
    plt.xlabel('$x$ or $y$ coordinate (units of $M$)')
    plt.ylabel('$z$ coordinate (units of $M$)')

def r_range(theta, num):
    r_plus = M + np.sqrt(M**2 - a**2)
    r_ergo = M + np.sqrt(M**2 - (a*np.cos(theta))**2)
    
    if r_plus*(1+1e-6) >= r_ergo*(1-1e-6):
        print(theta)
        raise ValueError
    return np.linspace(r_plus*(1+1e-6), r_ergo*(1-1e-6), num=num)
    

@np.vectorize
def best_alpha(r,theta):
    rhosq = r**2 + (a*np.cos(theta))**2
    
    g_tt = - (1 - 2 * M * r / rhosq)
    g_tp = - 2 * M * a * r * np.sin(theta)**2 / rhosq
    g_pp = (r**2 + a**2 + 2 * M * r * a**2 * np.sin(theta)**2 / rhosq) * np.sin(theta)**2
    
    def vsq(alpha): 
        return (
            alpha ** 2 * g_tt 
            + (1-alpha)**2 * g_pp 
            + 2 * alpha * (1-alpha) * g_tp
        )

    res = minimize_scalar(vsq, bracket=(0, 1), tol=1e-19)

    if res.success == True and res.fun < 0:
        return res.x
    else:
        print(res)
        print(r, theta)
        raise ValueError

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(kerr_killing)