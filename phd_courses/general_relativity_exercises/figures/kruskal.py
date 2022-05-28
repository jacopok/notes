import numpy as np
import matplotlib.pyplot as plt

def killing_kruskal():
    r_range = np.linspace(-3, 3)
    t_range = np.linspace(-3, 3)
    
    R, T = np.meshgrid(r_range, t_range)
    
    V = T + R
    U = T - R
    
    plt.quiver(R, T, 0, 2*R)
    
    plt.plot(r_range, r_range, c='black')
    plt.plot(r_range, -r_range, c='black')
    
    plt.plot(r_range, np.sqrt(r_range**2 + 1), c='black')
    plt.plot(r_range, -np.sqrt(r_range**2 + 1), c='black')
    
    plt.gca().set_aspect('equal')
    
    plt.xlabel('$R = (V - U) / 2$')
    plt.ylabel('$T = (V + U) / 2$')

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(killing_kruskal)