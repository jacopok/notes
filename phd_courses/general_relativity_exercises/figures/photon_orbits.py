import numpy as np
import matplotlib.pyplot as plt

def impact_parameter(Vmax):
    return 1/ (np.sqrt(Vmax)*2)

def effective_potential():
    rs = np.linspace(2, 8, num=1000)
    
    V = lambda r : 1 / 2 / r**2 - 1 / r**3
    
    plt.plot(rs, V(rs))
    
    for r in [3, 5]:
        less_rs = np.linspace(r, 8)
        plt.plot(less_rs, np.ones_like(less_rs) * V(r), label=impact_parameter(V(r)))

    plt.plot(rs, np.ones_like(rs)*1.2*V(3), label=impact_parameter(1.2*V(3)))
    
    plt.xlim(0, 8)
    plt.legend()
    

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(effective_potential)