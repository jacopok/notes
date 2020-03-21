M = .5 
omegak=.7 
V0 = np.linspace(0, 2*(M+omegak), num=1000) 
k = abs(((omegak-V0)**2 - M**2))**(1/2) 
plt.plot(V0, k) 
plt.xlabel('$V_0$') 
plt.ylabel('$\\abs{k_x^{\\prime}}$') 
plt.axvline(x=M, label='$M$', ls=':', c='red') 
plt.axvline(x=omegak, label='$\\omega_k$', ls=':', c='orange') 
plt.legend() 
plt.show(block=False) 
lsave('primed_momentum_potential')                                                          
