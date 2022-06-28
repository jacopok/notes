import numpy as np
import matplotlib.pyplot as plt
from matplotlib import transforms

def conformal_trajectory():
    lam = np.linspace(-5, 5, num=1000)
    a = 1
    t = np.sinh(a * lam) / a
    r = np.cosh(a * lam) / a
    
    u = t - r
    v = t + r
    
    U = np.arctan(u) 
    V = np.arctan(v)
    
    base = plt.gca().transData
    rot = transforms.Affine2D().rotate_deg(45)
    
    plt.plot(U, V, transform = rot + base)
    
    UV_range = np.linspace(-np.pi/2, np.pi/2, num=1000)
    plt.plot(UV_range, UV_range, transform = rot + base)
    plt.plot(UV_range, np.pi/2 * np.ones_like(UV_range), transform = rot + base)
    plt.plot(np.pi/2*np.ones_like(UV_range), UV_range, transform = rot + base)
    plt.plot(UV_range, -np.pi/2 * np.ones_like(UV_range), transform = rot + base)
    plt.plot(-np.pi/2*np.ones_like(UV_range), UV_range, transform = rot + base)

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(conformal_trajectory)