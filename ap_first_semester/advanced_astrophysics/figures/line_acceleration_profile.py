import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

rs = np.linspace(1, 3, num=400)

a = .6
b = .7

gs = rs ** (-2) * (1 - 1 / rs)**(a * (2 * b - 1))

plt.plot(rs, gs)
plt.xlabel("Normalized radius $r / R_*$")
plt.ylabel("Acceleration due to the line (arbitrary units) $g_L$")
plt.show(block=False)