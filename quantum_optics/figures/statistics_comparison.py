import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from scipy.special import factorial
plt.style.use(astropy_mpl_style)

def thermal(n, n_bar):
  return ((n_bar/(1+nbar)) ** n / (1 + n_bar))
  
def coherent(n, n_bar):
  return (np.exp(-n_bar) * n_bar ** n / factorial(n))
  
ns = np.arange(0, 40)
nbars = [6]
width = 0.35
alpha = .8

for nbar in nbars:
  plt.bar(ns, coherent(ns, nbar), label='coherent ' + str(nbar), alpha = alpha, width=width)
  plt.bar(ns+width, thermal(ns, nbar), label='thermal ' + str(nbar), alpha = alpha, width=width)
plt.legend()
plt.show()

def moment(b, v, n):
  m = np.average(b, weights=v)
  if (n == 1):
    return(m)
  return (np.sum((b - m)**n * v) / np.sum(v))

def analyze(dist, ns, nbar):
  bins = dist(ns, nbar)
  moments = {'mean': 1, 'variance': 2, 'skewness': 3, 'kurtosis': 4}
  for m, num in moments.items():
    print(m, ' = ', moment(ns, bins, num))

"""
For a Poissonian (coherent): mean=variance=skewness=nu 
and kurtosis = nu*(1+3*nu)

For the thermal distribution: nth moment is something like exp(n)
"""
