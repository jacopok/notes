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

from basic_units import cos,sin, radians

theta_r = np.linspace(0, 2 * np.pi, num=1000)
theta = np.array([x * radians for x in theta_r])

plt.plot(theta, (1 - cos(theta)) / 2, label='Scale factor $a$')
plt.plot(theta, (theta_r - sin(theta)) / 2, label='Time $t$')
plt.xlabel('Angle parameter $\\theta$')
plt.ylabel('Normalized scale factor $a / \\widetilde{a}_0$ or time $t / \\widetilde{t}_0$')
plt.legend()
plt.savefig('positive_curvature_a.pdf', format = 'pdf')
# plt.show()