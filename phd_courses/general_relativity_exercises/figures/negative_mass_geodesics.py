import numpy as np
import matplotlib.pyplot as plt
cmap = plt.get_cmap('plasma')


def negative_mass_geodesics():
    rs = np.logspace(-1, 1, num=1000)
    
    M = -1
    
    V_null = lambda r : 1 / 2 / r**2 - M / r**3
    V_mass_l = lambda r : - M / r + 1/2/r**2 - M / r**3
    V_mass_nol = lambda r : - M / r
    
    plt.loglog(rs, V_null(rs), c=cmap(.1), label='Null, nonzero $l$')
    plt.loglog(rs, V_mass_l(rs), c=cmap(.5), label='Timelike, nonzero $l$')
    plt.loglog(rs, V_mass_nol(rs), c=cmap(.9), label='Timelike, zero $l$')
    

    plt.xlabel('Radius [units of $M$]')
    plt.ylabel('Potential [dimensionless]')
    plt.yticks(None)
    plt.legend()
    

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(negative_mass_geodesics)