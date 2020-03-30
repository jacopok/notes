import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)


ax = plt.subplot(projection='polar')
ax.set_theta_zero_location("N")

thetas = np.linspace(0, 2 * np.pi, num=1000)

plus_polarization_amplitude = (1 + np.cos(thetas)**2)/2
cross_polarization_amplitude = np.cos(thetas)

ax.plot(thetas, ((plus_polarization_amplitude)**2 + (cross_polarization_amplitude)**2) * abs(np.sin(thetas))) 
plt.show(block=False)