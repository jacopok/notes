import numpy as np
import matplotlib.pyplot as plt


def q_function(beta_r, beta_i, alpha, r, theta):
    beta = beta_r + 1j * beta_i

    return (
        1/(np.pi * np.cosh(r))
        * np.exp(
            -(abs(alpha)**2 + abs(beta)**2)
            + (
                np.conj(beta) * alpha +
                np.conj(alpha) * beta
            ) / np.cosh(r)
            - 1/2 * (
                np.exp(1j * theta) * (
                    np.conj(beta)**2 - np.conj(alpha)**2
                )
                + np.exp(-1j*theta) * (
                    beta**2 - alpha**2
                )
            ) * np.tanh(r)
        )
    ).real


def husini():
    num = 1000
    beta_r_range = np.linspace(-4, 4, num=num)
    beta_i_range = np.linspace(-4, 4, num=num)
    
    beta_r, beta_i = np.meshgrid(beta_r_range, beta_i_range)
    
    params = {
        (0, 0): {'alpha': 0, 'r': 0, 'theta': 0},
        (0, 1): {'alpha': 0, 'r': 1, 'theta': 0},
        (1, 0): {'alpha': 2, 'r': 0, 'theta': np.pi/2},
        (1, 1): {'alpha': 1+1j, 'r': 1, 'theta': np.pi/2},
    }
    
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, squeeze=False)
    
    for idx, ax in np.ndenumerate(axs):
    
        q = q_function(beta_r, beta_i, **params[idx])
        color_mappable = ax.contourf(beta_r, beta_i, q, levels=100)
        ax.set_title(
            f'$\\alpha=${params[idx]["alpha"]:.1f}, '
            f'$\\xi=${params[idx]["r"]*np.exp(1j*params[idx]["theta"]):.1f}'
        )
        
    cbar = fig.colorbar(color_mappable, ax=axs[:,:])
    cbar.ax.set_ylabel('Husini $Q$ function')

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel('Real part of $\\beta$')
    plt.ylabel('Imaginary part of $\\beta$')

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(husini)
