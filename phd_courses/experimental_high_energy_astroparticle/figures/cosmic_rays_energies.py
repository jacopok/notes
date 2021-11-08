import numpy as np
import matplotlib.pyplot as plt
          
energies = np.logspace(8, 20)
THR = 1e15
SCALING = min(energies) ** (2.7)

@np.vectorize
def flux(E):
    if E < THR:
        return E**(-2.7) * (THR)**(2.7-3) * SCALING
    else:
        return E**(-3) * SCALING

def cosmic_rays_energies():
    plt.loglog(energies, flux(energies))
    plt.grid('on')
    plt.xlabel('Energy [eV]')
    plt.ylabel('Flux [arbitrary scaling]')
