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

a = (1 - cos(theta)) / 2
t = (theta_r - sin(theta)) / 2

plt.plot(theta, a, label='Scale factor $a$')
plt.plot(theta, t, label='Time $t$')
plt.xlabel('Angle parameter $\\theta$')
plt.ylabel('Normalized scale factor $a / \\widetilde{a}_0$ or time $t / \\widetilde{t}_0$')
plt.legend()
plt.savefig('positive_curvature_a.pdf', format = 'pdf')
# plt.show()
plt.close()

plt.plot(t, a)
plt.xlabel("Time $t / \\widetilde{t}_0$")
plt.ylabel("Scale factor $a / \\widetilde{a}_0$")
plt.savefig('positive_curvature_a_vs_t.pdf', format = 'pdf')