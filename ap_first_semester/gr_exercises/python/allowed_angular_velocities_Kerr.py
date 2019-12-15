import numpy as np
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

def omega1(r):
    omega1 = (2 - r) / (2 - r + r ** 2)
    return(omega1)

def omega2(r):
    omega2 = 1 / (1 + r)
    return (omega2)

rs = np.linspace(1, 6, num=200)
plt.plot(rs, omega1(rs), label="$\\Omega_{-}$")
plt.plot(rs, omega2(rs), label="$\\Omega_{+}$")
plt.xlabel("$R = r/GM$")
plt.ylabel("$O = \\Omega GM$")
plt.legend()
plt.savefig(fname='../figures/allowed_velocities.pdf', format = 'pdf')