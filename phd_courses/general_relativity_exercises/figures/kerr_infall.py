import numpy as np
import matplotlib.pyplot as plt

r = np.linspace(2, 10, num=1000)

e = .001
l = 0
a = .5
M = 1

rhosq = r**2 + a**2
delta = r**2 - 2*M*r + a**2

g_tt = (a**2 - delta) / rhosq
g_tp = a * (delta / rhosq - 1)
g_pp = rhosq - a**2 * delta / rhosq
g_rr = rhosq / delta

tdot = e / (g_tp**2 / g_pp - g_tt)
pdot = - g_tp / g_tt * tdot

rdotsq = ( -1 + e * tdot) / g_rr

plt.plot(r, tdot, label=r'$\dot{t}$')
plt.plot(r, pdot, label=r'$\dot{\varphi}$')
plt.plot(r, np.sqrt(rdotsq), label=r'$|\dot{r}|$')
plt.legend()
plt.show()