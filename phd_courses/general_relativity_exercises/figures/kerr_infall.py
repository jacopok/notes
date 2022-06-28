import numpy as np
import matplotlib.pyplot as plt

r = np.linspace(2, 10, num=1000)

e = .1
l = 0
a = .5
M = 1

theta = np.pi / 2

rhosq = r**2 + (a*np.cos(theta))**2
delta = r**2 - 2*M*r + a**2

g_tt = - (1 - 2 * M * r / rhosq)
g_tp = - 2 * M * a * r * np.sin(theta)**2 / rhosq
g_pp = (r**2 + a**2 + 2 * M * r * a**2 * np.sin(theta)**2 / rhosq) * np.sin(theta)**2
g_rr = rhosq / delta

tdot = e / (g_tp**2 / g_pp - g_tt)
pdot = - g_tp / g_tt * tdot

rdotsq = ( -1 + e * tdot) / g_rr

plt.plot(r, tdot, label=r'$\dot{t}$')
plt.plot(r, pdot, label=r'$\dot{\varphi}$')
plt.plot(r, rdotsq, label=r'$|\dot{r}|^2$')
plt.legend()
plt.show()