import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import matplotlib.animation as animation
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('text.latex')
        

thetas = np.linspace(0, 2 * np.pi, num=1000)

LL = (1 + np.cos(thetas))**2
LR = (1 - np.cos(thetas))**2

fig = plt.figure()
ax = fig.add_subplot(projection='polar')
ax.set_theta_zero_location("N")

ax.plot(thetas, LL, ls =':')
ax.plot(thetas, LR, ls = '--')

# ax.legend()
ax.get_yaxis().set_visible(False)
plt.tight_layout()
fig.savefig('angular_distribution.png', format = 'png', dpi=300)