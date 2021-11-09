import matplotlib.pyplot as plt
import numpy as np
from functools import partial

GAMMA_LIST = [1, 2, 10]
BETA_PRODUCT = .5
NUM = 200

def color_num(gamma):
    return ((np.log(gamma) - np.log(min(GAMMA_LIST))) 
    / (np.log(max(GAMMA_LIST)) - np.log(min(GAMMA_LIST))))

def beta(gamma):
    return np.sqrt(1 - 1 / gamma)

def boost_angle_com_to_lab(theta, gamma, beta_product):
    return np.arctan(
        np.sin(theta)
        / gamma
        / (
            np.cos(theta)
            + beta(gamma) / beta_product
        )
    ) % (np.pi)

def angular_distribution_shift():

    angles_com = np.linspace(0, np.pi, num=NUM)
    weights_com = np.sin(angles_com)
    
    colormap = plt.get_cmap('viridis')
    def color(gamma):
        return colormap(color_num(gamma))

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    
    for gamma in GAMMA_LIST:
        angles_lab = boost_angle_com_to_lab(angles_com, gamma=gamma, beta_product=BETA_PRODUCT)
        ax.plot(angles_lab, weights_com, label=f'$\gamma=${gamma}', color=color(gamma))
        ax.plot(2 * np.pi - angles_lab, weights_com, color=color(gamma))
    ax.set_theta_zero_location("N")
    ax.legend()
    ax.set_yticklabels([])

if __name__ == "__main__":

    from make_all_figures import plot_and_save
    plot_and_save(angular_distribution_shift)
