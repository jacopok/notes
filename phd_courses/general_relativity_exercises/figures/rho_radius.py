import numpy as np
import matplotlib.pyplot as plt

def rho_radius():

    rho = np.linspace(1/4, 3)

    r = rho * (1 + 1 / 2 / rho)**2

    print(min(r))
    print(min(rho))

    plt.plot(r, rho)


    

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(rho_radius)