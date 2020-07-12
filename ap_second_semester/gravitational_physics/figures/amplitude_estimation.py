import numpy as np
from astropy.constants import codata2018 as ac
import astropy.units as u

m1 = 36 * u.Msun
m2 = 29 * u.Msun
m = (m1 * m2 ) / (m1 + m2)
M = m1 + m2             

f_gw = 35 * u.Hz
omega_s = (2 * np.pi * f_gw) / 2

r = 410 * u.Mpc
A = (4 * ac.G**(5/3) * (omega_s)**(2/3) * m * M**(2/3) / r / ac.c**4).to(1)
