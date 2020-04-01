import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
plt.xkcd(randomness=1)

R1 = np.linspace(0, 1)
R2 = np.linspace(1, 6)

v_ob1 = R1
v_ob2 = np.ones_like(R2)
v_th2 = R2 ** (-1 / 2)

plt.plot(R1, v_ob1, c='black')
plt.plot(R2, v_ob2, c='black', ls = '-', label='Observed')
plt.plot(R2, v_th2, c='black', ls='--', label='Predicted')

plt.legend()
# plt.axis('off')
plt.ylabel('$v(R)$')
plt.xlabel('$R$')

plt.show()