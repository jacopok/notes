import astropy.units as u
from astropy.constants import codata2018 as ac

@u.quantity_input(M='mass', D='length', v='speed')
def  estimate_h(M, D, v) -> u.dimensionless_unscaled:
    h = ac.G * M / D / ac.c**2 * (v / ac.c)**2
    # the units system takes care of all the unit conversion for us
    # since we specified that we want the result to be a pure number
    return(h)

print(f'Car crash: {estimate_h (1e3*u.kg, 10*u.m, 100* u.km/u.hr):.0e}')
print(f'Supernova: {estimate_h (1*u.Msun , 2*u.kpc , .2 * ac.c):.0e}')
print(f'Binary BH: {estimate_h (50*u.Msun , 400*u.Mpc , .1 * ac.c):.0e}')