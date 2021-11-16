import numpy as np
import matplotlib.pyplot as plt

NUM = int(1e4)

def pseudo_rapidity(theta: np.ndarray):
    return - np.log(np.arctan(theta / 2))

def rapidity():
    thetas = np.linspace(0, np.pi, num=NUM)
    plt.plot(thetas, pseudo_rapidity(thetas))
    plt.grid('on')

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(rapidity)