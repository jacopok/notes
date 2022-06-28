import numpy as np
import matplotlib.pyplot as plt
from matplotlib import transforms

cmap = plt.get_cmap('plasma')

def conformal_trajectory():
    lam = np.linspace(-20, 20, num=1000)
    
    base = plt.gca().transData
    rot = transforms.Affine2D().rotate_deg(45)
    
    for a in [1/2, 1, 2]:
        t = np.sinh(a * lam) / a
        r = np.cosh(a * lam) / a
        
        u = t - r
        v = t + r
        
        U = np.arctan(u)
        V = np.arctan(v)
        
        plt.plot(V, U, transform = rot + base, label = f'a={a}', color=cmap(a/3))
    
    ph = np.pi/2
    
    UV_range = np.linspace(-ph, ph, num=1000)
    kwargs = {
        'transform': rot + base,
        'color': 'black'
    }
    
    kwargs_curve = {
        'transform': rot + base,
        'color': 'green'
    }
    
    plt.plot(UV_range, UV_range, **kwargs, ls='--')
    plt.plot(ph*np.ones_like(UV_range), UV_range, **kwargs)
    plt.plot(UV_range, -ph * np.ones_like(UV_range), **kwargs)

    t = np.linspace(- 100, 0, num=1000)
    r = np.zeros_like(t)
    u = t - r
    v = t + r    
    U = np.arctan(u)
    V = np.arctan(v)    
    plt.plot(V, U, **kwargs_curve)

    t = np.linspace(0, 1, num=200)
    r = np.sqrt(1 + t**2) - 1
    u = t - r
    v = t + r    
    U = np.arctan(u)
    V = np.arctan(v)
    plt.plot(V, U, **kwargs_curve, ls=':')

    t = np.linspace(1, 100, num=200)
    r = t / np.sqrt(2) - 1/ np.sqrt(2) + np.sqrt(2) - 1
    u = t - r
    v = t + r    
    U = np.arctan(u)
    V = np.arctan(v)    
    plt.plot(V, U, **kwargs_curve)


    plt.text(ph+.1, ph-.05, '$\iota^+$', **kwargs)
    plt.text(ph+.1, -.05, r'$I^+$', **kwargs)
    plt.text(-ph+.1, -ph-.1, '$\iota^-$', **kwargs)
    plt.text(ph+.1, -ph-.05, '$\iota^0$', **kwargs)
    # plt.plot(UV_range, ph * np.ones_like(UV_range), **kwargs)
    # plt.plot(-ph*np.ones_like(UV_range), UV_range, **kwargs)
    plt.gca().set_aspect('equal')
    plt.legend()

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(conformal_trajectory)