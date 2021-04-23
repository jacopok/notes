from astropy.constants import codata2018 as ac
import astropy.units as u
import numpy as np

M = 1.4 * u.Msun
h = 1.6 * u.m
R = 12 * u.km

@u.quantity_input
def acc_newt(r) -> u.m / u.s**2:
    return ac.G * M / r**2 

@u.quantity_input
def acc_gr(r) -> u.m / u.s**2:
    return ac.G * M / r**2 / np.sqrt(1 - 2 * ac.G * M / ac.c**2 / r)

print('Newtonian, approximated')
print((2 * ac.G * M * h / R**3).to(u.m / u.s**2))

print('Newtonian, exact')
print(acc_newt(R) - acc_newt(R+h))

print('GR, exact')
print(acc_gr(R) - acc_gr(R + h))

print('Non-free-fall ratio')
print((acc_gr(R) - acc_gr(R + h)) / (acc_newt(R) - acc_newt(R+h)))