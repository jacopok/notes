import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('text.latex', preamble=r'''\usepackage{amsmath}
          \usepackage{physics}
          \usepackage{siunitx}
          ''')
rc('figure', dpi=150)
          
energies = np.logspace(8, 20)

@np.vectorize
def flux(E):
    if E < 1e15:
        return E**(-2.7)
    else:
        return E**(-3)

if __name__ == "__main__":
    plt.loglog(energies, flux(energies))
    plt.savefig(str(__file__).split(sep='.')[0] + '.pdf')