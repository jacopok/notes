import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text.latex', preamble=r'''\usepackage{amsmath}
\usepackage{physics}
\usepackage{siunitx}
''')
R = np.linspace(1, 1.5, num=200)
f = lambda x: np.arcsin(np.sqrt(27/4/x**2 * (1 - 1/x))) / (np.pi/2)
fig = plt.figure(1)
plt.plot(R, f(R))
plt.grid()
# plt.show()
plt.xlabel('$R = r / 2GM$')
plt.ylabel('$\psi_{\\text{critical}} / (\\pi/2)$')
plt.savefig('../figures/critical_psi.pdf', format = 'pdf')
plt.close()