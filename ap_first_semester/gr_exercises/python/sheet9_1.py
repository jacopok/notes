import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('text.latex', preamble=r'''\usepackage{amsmath}
          \usepackage{physics}
          \usepackage{siunitx}
          ''')

import numpy as np

def V(e, r):
    t1 = 3/(8*r**2) - 1/(8*r**3)
    t2 = - 1/(2*r) + 1/(8*r**2)
    return(e**2 * t1 + t2)

x= np.linspace(1/3, 3, num=300)

for e in np.arange(1, 3, .5):
    plt.plot(x, V(e, x) - (e**2-1)/2, label=f"$e =$ {e}")
plt.legend()
plt.xlabel("$R$")
plt.ylabel("$ V _{\\text{eff}} (R, e) - (e^2-1)/2$")

plt.savefig('../figures/potential_barrier.pdf', format='pdf')
plt.close()