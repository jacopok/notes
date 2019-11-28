import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from matplotlib import rc
from scipy.optimize import fsolve
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text.latex', preamble=r'''\usepackage{amsmath}
\usepackage{physics}
\usepackage{siunitx}
''')

plt.close()

R = np.linspace(0, 3, num=200)
f = lambda x: (x-1) * np.exp(x/2) 
# fig = plt.figure(1)
# plt.plot(R, f(R))
# plt.grid()
# plt.show()
# plt.close()

@np.vectorize
def r(U, V, hint=None):
    """
    a numerical solution to 
    U**2 - V**2 = (r-1) * exp(r/2)
    """

    LHS = U ** 2 - V ** 2
    if (LHS < 0):
        return (0)
    function = lambda x: (x-1)*np.exp(x/2) - LHS
    if(not hint):
        hint = max(np.log(LHS), 0)
    
    return fsolve(function, x0=hint)
        
Us = np.linspace(-10, 10, num=200)
for V0 in range(4):
    rs = r(Us, V0)
    if (V0 == 0):
        lab = f"$V_0 = $ {V0}"
    else:
        lab= f"$V_0 = \pm$ {V0}"
    plt.plot(rs, Us, label=lab)
plt.xlabel("$r(U^2 - V^2) / 2GM$")
plt.ylabel("$U$")
plt.legend()
plt.show()