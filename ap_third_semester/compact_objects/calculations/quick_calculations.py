#%%

import numpy as np

print(1 - np.sqrt(8 / 9))
print(1 / 12)

#%%

theta = np.linspace(0, 2* np.pi)
Cs = np.linspace(0, 2, num=10)

ax = plt.subplot(projection='polar')
ax.set_theta_zero_location('N')

for C in Cs:
    ax.plot(theta, C * np.sin(theta)**2)

# %%
