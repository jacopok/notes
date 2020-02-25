import scipy.constants as sc
from astropy.constants import codata2018 as ac
from astropy.constants import iau2015 as aa
import astropy.units as u
from astropy.cosmology import Planck15 as cosmo
import astropy.uncertainty as aun
import pandas as pd
from functools import partial

names = ['Hydrogen', 'Helium', 'Carbon', 'Neon', 'Oxygen', 'Silicon']
temperatures = [3, 14, 53, 110, 160, 270] * u.keV
L_fractions = [0, 1.7e-5, 1, 1.8e3, 2.1e4, 9.2e5] * u.dimensionless_unscaled
L_gammas = [2.1, 6, 8.6, 9.6, 9.6, 9.6] * u.Unit(u.L_sun * 1e4)
rhos = [5.9, 1.3e3, 1.7e5, 1.6e7, 9.7e7, 2.3e8] * u.g / u.cm**3
durations = [1.2e7, 1.3e6, 6.3e3, 7, 1.7, 6 / 365.25] * u.yr

start_el = ['H', 'He', 'C', 'Ne', 'O', 'Si']
end_el = ['He', 'C, O', 'Ne, Mg', 'O, Mg', 'Si', 'Fe, Ni']
processes = []

for s, e in zip(start_el, end_el):
  processes.append(s + ' $\\rightarrow$ ' + e)

cols = [names, processes, temperatures, rhos, L_gammas, L_fractions, durations]

pd_cols = []
for x in cols:
  pd_cols.append(pd.DataFrame(x))

df = pd.concat(pd_cols, axis=1)

def si_format(x, digits=0):
  if(isinstance(x, u.Quantity)):
    return(f'\\num{{{x.value:.{digits}e}}}')
  else:
    return (f'\\num{{{x:.{digits}e}}}')

def float_format(x, digits=2):
  if(isinstance(x, u.Quantity)):
    return(f'\\num{{{x.value:.{digits}f}}}')
  else:
    return(f'\\num{{{x:.{digits}f}}}')

int_format = partial(float_format, digits=0)
one_digit = partial(float_format, digits=1)
one_digit_exp = partial(si_format, digits=1)


ident = lambda x : x

formatters = [ident, ident, int_format, one_digit_exp, one_digit, one_digit_exp, one_digit_exp]

# print(df.to_latex(formatters=formatters, escape=False, index=False))