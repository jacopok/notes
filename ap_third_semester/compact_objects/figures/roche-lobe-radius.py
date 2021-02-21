import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('text.latex', preamble=r'''\usepackage{amsmath}
          \usepackage{physics}
          \usepackage{siunitx}
          ''')
THR = .5
WIDTH = 0

# def weight(q):
#     if WIDTH>0:
#         offset = 1/2 - THR / WIDTH
#         return (np.piecewise(q, 
#             condlist=[
#                 q < THR - WIDTH / 2,
#                 q > THR - WIDTH / 2 and q < THR + WIDTH / 2 ,
#                 q > THR + WIDTH / 2,
#             ],
#             funclist=[
#                 0,
#                 lambda x: x / WIDTH + offset,
#                 1
#             ]
#         ))
#     else:
#         return (np.piecewise(q,
#             condlist=[q < THR, q >= THR],
#             funclist=[0, 1]
#         ))

def f1(q):
    return (.46224 * (q / (1 + q))**(1 / 3))
    
def f2(q):
    return (.38 + .2 * np.log10(q))


def f(q):
    if q < 0.5:
        return (f1(q))
    else:
        return(f2(q))

f = np.vectorize(f, signature='()->()')

qs = np.linspace(0, 8, num=1000)

f_q = f(qs)

def a(q):
    return((1+q)**4 / q**2)

a_q = a(qs)

plt.plot(qs, np.abs(np.gradient(f_q, qs) / f_q), label='$\\abs{\\Delta \\log f}$')
plt.plot(qs, np.abs(np.gradient(a_q, qs) / a_q), label='$\\abs{\\Delta \\log a}$')
plt.plot(qs, np.gradient(a_q, qs) / a_q + np.gradient(f_q, qs) / f_q, label='$\\Delta \\log a + \\Delta \\log f$', ls='--')
plt.axvline(1, label='$q = 1$', ls=':', c='black')

# plt.plot(qs, f(qs))
# plt.xlabel('$q = M_2 / M_1$')
# plt.ylabel('$R_{\\text{{lobe}}} / a$')
# plt.savefig('roche-lobe-radius.pdf', format = 'pdf')

plt.xlabel('$q = M_2 / M_1$')
plt.ylabel('relative variation')
plt.legend()
plt.yscale('log')
plt.savefig('roche-lobe-relative-corrections.pdf')
plt.show()