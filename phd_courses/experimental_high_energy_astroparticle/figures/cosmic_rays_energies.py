import numpy as np
import matplotlib.pyplot as plt
          
energies = np.logspace(-1, 11)
THR = 1e6
SCALING = 1.8e4

# Formula from PDG 2020, page 219

@np.vectorize
def flux(E):
    if E < THR:
        return E**(-2.7) * (THR)**(2.7-3) * SCALING
    else:
        return E**(-3) * SCALING

def cosmic_rays_energies():
    plt.loglog(energies, flux(energies))
    plt.grid('on')
    plt.xlabel('Energy [GeV]')
    plt.ylabel('Flux [$\SI{}{\meter^{-2}\ \second^{-1}\ \steradian^{-1}\ \GeV^{-1}}$]')

if __name__ == "__main__":

    from make_all_figures import plot_and_save
    plot_and_save(cosmic_rays_energies)
