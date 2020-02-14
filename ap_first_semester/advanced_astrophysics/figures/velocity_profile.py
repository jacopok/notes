import numpy as np
import matplotlib.pyplot as plt
plt.style.use('default')
from scipy.optimize import fsolve
import astropy.units as u
import astropy.constants as ac

def Mach(r, r0, rc, a, v0):
  RHS = ((v0 / a / (r / r0)** 2 * np.exp(2 * rc * (1 / r0 - 1 / r)))).to('').value
  eq = lambda x: x * np.exp(-x ** 2 / 2) - RHS
  if (r > rc):
    guess = 1.5
  else:
    guess = 0.5
  M = fsolve(eq, guess)
  return (M[0])

def density(r, r0, rc, a, v0):
  """
  Returns rho / rho0
  """
  RHS = (np.exp(-2 * rc * (1 / r0 - 1 / r))).value
  exp_arg = (1/2 * (v0 * r0**2 / a / r**2)**2).to('').value
  eq = lambda x: x * np.exp(x ** (-2) * exp_arg) - RHS
  if (r < rc):
    guess = 1e-3
  else:
    guess = (atmosphere_density(r, r0, rc) / (r / rc)**2 /2).to('')
  print(guess)
  rho = fsolve(eq, guess)
  return (rho[0])
  
def atmosphere_density(r, r0, rc):
  return(np.exp(-2 * rc *(1/r0 - 1/r)))
  
r0 = 1 * u.R_sun
rs = np.logspace(0, 3) * r0
R = ac.k_B / u.u
mu = 0.62
T = 1e6 * u.K
a = np.sqrt(R * T / mu)
rc = ac.GM_sun / 2 / a ** 2
v0 = a * (rc / r0)**2 * np.exp(-2*rc/r0 + 3/2)

Ms = []
for r in rs:
  Ms.append(Mach(r, r0, rc, a, v0))

rhos = []
atm_rhos = []
for r in rs:
  rhos.append(density(r, r0, rc, a, v0))
  atm_rhos.append(atmosphere_density(r, r0, rc))

fig, axs = plt.subplots(2, 1, sharex=True)
axs[0].loglog(rs, Ms)
axs[0].set_ylabel("Mach number: $v / a = M$")
axs[1].loglog(rs, rhos)
axs[1].loglog(rs, atm_rhos, linestyle = ':', linewidth=2., label="Atmosphere", c='black')
axs[1].set_xlabel("Radius: $r/r_0$")
axs[1].set_ylabel("Density: $\\rho / \\rho_0$")
for ax in axs:
  ax.axvline(rc / r0, label = "Critical radius", linestyle="--", c="purple")
  ax.grid('on')
plt.legend()
plt.show(block=False)
