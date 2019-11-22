import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style

R = np.linspace(1, 5, num=200)

f = lambda x: x ** (-2) - x ** (-3)

fig = plt.figure(1)
plt.plot(R, f(R), label = 'Function 1/R^2 - 1/R^3')
plt.xlabel('R = r/r_s')
plt.legend()
plt.grid()
plt.savefig(fname="photon_effective_potential.pdf", format="pdf")
