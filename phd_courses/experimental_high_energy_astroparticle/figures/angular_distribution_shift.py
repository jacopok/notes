import matplotlib.pyplot as plt
import numpy as np
from functools import partial
from scipy import stats

GAMMA_LIST = np.logspace(0, .8, num=4)
BETA_PRODUCT = 1
NUM = 20_000_000
BIN_NUMBER = 200

LOG_FUDGE_FACTOR = .2

def color_num(gamma):
    return ((np.log(gamma) - np.log(min(GAMMA_LIST))) 
    / (
        np.log(max(GAMMA_LIST)) 
      - np.log(min(GAMMA_LIST))
      + LOG_FUDGE_FACTOR)
    )

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

class SineRandomVariable(stats.rv_continuous):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.a = 0.
        self.b = np.pi

    def _pdf(self, theta):
        return (np.sin(theta) / 2)
    
    def _ppf(self, x):
        return(np.arccos(1 - 2*x))


def angular_distribution_shift():

    angles_com = SineRandomVariable().rvs(size=NUM)
    
    colormap = plt.get_cmap('plasma')
    def color(gamma):
        return colormap(color_num(gamma))

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    
    for gamma in GAMMA_LIST:
        angles_lab = boost_angle_com_to_lab(angles_com, gamma=gamma, beta_product=BETA_PRODUCT)
        
        vals, bins_with_edges = np.histogram(angles_lab, bins=BIN_NUMBER)
        bins = (
            bins_with_edges[1:]+ 
            bins_with_edges[:-1]
            ) / 2
        vals = vals / max(vals)
        
        ax.plot(bins, vals, label=f'$\gamma=$ {gamma:.1f}', color=color(gamma))
        # ax.plot(2 * np.pi - bins, vals, color=color(gamma))
        
        ax.axvline(1 / gamma, color=color(gamma), ls=':')
        # ax.axvline(-1 / gamma, color=color(gamma), ls=':')
        
        # fraction_below_gamma_inv = (
        #     sum(1 for ang in angles_lab if ang < 1/ gamma)
        #     / len(angles_lab)
        # )
        
        # print(f'Fraction sent below $1/\\gamma$: {100*fraction_below_gamma_inv:.1f}%')
        
    # ax.set_theta_zero_location("N")
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    # ax.set_xlabel('$\\theta$')
    # ax.set_ylabel('$p(\\theta)$')
    ax.legend()
    ax.set_yticklabels([])
    plt.tight_layout()
    

if __name__ == "__main__":

    from make_all_figures import plot_and_save
    plot_and_save(angular_distribution_shift)
