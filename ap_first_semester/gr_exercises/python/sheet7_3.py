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
R = np.linspace(1.5, 5, num=200)
f = lambda x: np.sqrt(1 + 1/(2*x-3))/np.sqrt(1 - 1/x)
fg = lambda x: 1 / np.sqrt(1 - 1 / x) - 1
fd = lambda x: np.sqrt(1 + 1/(2*x-3))- 1
fig = plt.figure(1)
# plt.plot(R, f(R))
# plt.plot(R, f2(R))
plt.plot(R, fg(R)/fd(R))
plt.grid()
# plt.close()
