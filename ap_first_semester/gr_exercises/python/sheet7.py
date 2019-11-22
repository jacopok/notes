import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style

R = np.linspace(1, 5, num=200)

f = lambda x: x ** (-2) - x ** (-3)

fig = plt.figure(1)
plt.plot(R, f(R))
plt.xlabel('R = r/r_s')
plt.grid()
plt.savefig(fname="../figures/photon_effective_potential.pdf", format="pdf")
