import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import matplotlib.animation as animation
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('text.latex', preamble=r'''\usepackage{amsmath}
          \usepackage{physics}
          \usepackage{siunitx}
          ''')

NMAX = 200

thetas = np.linspace(0, 2 * np.pi, num=1000)

plus_polarization_amplitude = (1 + np.cos(thetas)**2)/2
cross_polarization_amplitude = np.cos(thetas)

fig = plt.figure()
ax = fig.add_subplot(projection='polar')
ax.set_theta_zero_location("N")
# line, = ax.plot([], [])

def init():
  ax.set_ylim((0,1.1))
  line.set_data([], [])
  return (line,)

def total_amplitude(time, Nmax=NMAX):
  t = np.pi * time / Nmax
  amplitude = (((plus_polarization_amplitude * np.sin(t))** 2 + (cross_polarization_amplitude * np.cos(t))** 2))
  
  line.set_data(thetas, amplitude)

def averaged_amplitude():
  amplitude = (((plus_polarization_amplitude)** 2 + (cross_polarization_amplitude)** 2))
  return(amplitude)

# anim = animation.FuncAnimation(fig, total_amplitude, range(NMAX), init_func=init, interval=50)

# anim.save('angular_spectrum_no_sin.gif', writer='imagemagick', fps=60, dpi=200)

ax.plot(thetas, averaged_amplitude())
ax.set_title('$r(\\theta) \\propto \\dv{E}{\\Omega}$')

plt.show(block=False)
