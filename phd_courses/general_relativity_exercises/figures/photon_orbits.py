import numpy as np
import matplotlib.pyplot as plt
cmap = plt.get_cmap('plasma')

def impact_parameter(Vmax):
    return 1/ (np.sqrt(Vmax*2))

# Vmax = e**2 / l**2 / 2 = 1/2 / b**2
# b**2 = 1 / 2 / Vmax

def effective_potential():
    rs = np.linspace(2, 8, num=1000)
    
    V = lambda r : 1 / 2 / r**2 - 1 / r**3
    
    plt.plot(rs, V(rs), c='black')
    
    plt.plot(rs, np.ones_like(rs)*1.2*V(3), 
             label=f'$b=${impact_parameter(1.2*V(3)):.2f}$M$',
             c=cmap(.1))
    for r in [3, 5]:
        less_rs = np.linspace(r, 8)
        plt.plot(less_rs, np.ones_like(less_rs) * V(r), 
                 label=f'$b=${impact_parameter(V(r)):.2f}$M$',
                 c=cmap(r/6))

    plt.xlim(0, 8)
    plt.xlabel('Radius [units of $M$]')
    plt.ylabel('Potential divided by $l^2$ [units of $M^{-2}$]')
    plt.yticks(None)
    plt.legend()
    

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(effective_potential)