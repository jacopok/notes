import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import matplotlib.animation as animation

NMAX = 200

thetas = np.linspace(0, 2 * np.pi, num=1000)

plus_polarization_amplitude = (1 + np.cos(thetas)**2)/2
cross_polarization_amplitude = np.cos(thetas)

fig = plt.figure()
ax = fig.add_subplot(projection='polar')
ax.set_theta_zero_location("N")
line, = ax.plot([], [])

def init():
  ax.set_ylim((0,1.1))
  line.set_data([], [])
  return (line,)

def total_amplitude(time, Nmax=NMAX):
  t = np.pi * time / Nmax
  amplitude = (((plus_polarization_amplitude * np.sin(t))** 2 + (cross_polarization_amplitude * np.cos(t))** 2))
  
  line.set_data(thetas, amplitude)

anim = animation.FuncAnimation(fig, total_amplitude, range(NMAX), init_func=init, interval=50)

anim.save('angular_spectrum_no_sin.gif', writer='imagemagick', fps=60, dpi=200)

# plt.show(block=False)
