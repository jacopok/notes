import numpy as np

import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from mpl_toolkits import mplot3d

x = np.linspace(-2, 2, num=20)
y = np.linspace(-2, 2, num=20)
X, Y = np.meshgrid(x, y)

def hyper(x, y):
  # return(np.sqrt(5 + x ** 2 - y ** 2))
  return(x ** 2 - y ** 2)

"""
Metric is:
g_μν =
[[1 + 4x**2, -4xy],
[-4xy, 1+4y**2]]

"""

Z = hyper(X, Y)

def plot_line(x0, y0, x1, y1):
  x_line = np.linspace(x0, x1)
  y_line = np.linspace(y0, y1)
  ax.plot(x_line, y_line, hyper(x_line, y_line), c='b')


fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='black', linewidth=.5, alpha=.6)
plot_line(-1.5, -1.5, 0.5, 0)
plt.show()