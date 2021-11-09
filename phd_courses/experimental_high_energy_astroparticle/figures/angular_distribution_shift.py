import matplotlib.pyplot as plt
import numpy as np

GAMMA_LIST = [1, 2, 10]

def boost_angle_com_to_lab(theta, gamma):

    pass    

def angular_distribution_shift():

    angles_com = np.linspace(0, np.pi)
    weights_com = np.sin(angles_com)
    
    for gamma in GAMMA_LIST:
        angles_lab = boost_angle_com_to_lab(angles_com)
        plt.plot(angles_lab, weights_com)

if __name__ == "__main__":

    from make_all_figures import plot_and_save
    plot_and_save(angular_distribution_shift)
