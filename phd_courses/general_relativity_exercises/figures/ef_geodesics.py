import numpy as np
import matplotlib.pyplot as plt

CMAP = plt.get_cmap('plasma')

r = np.linspace(0, 4, num=1000)

r0 = 4
e = np.sqrt(1 - 2 / r0)

c1 = CMAP(.8)
c2 = CMAP(.2)

rdot = - np.sqrt(e**2 + 2 / r - 1)
vdot = 1 / (e - rdot)
integral_coordinate = + ((1 + np.sqrt(-1 + 4/r)) * r) + 4 * np.arctan((np.sqrt(-1 + 4/r) * (-2 + r))/(-4 + r)) + 2 * np.log(2 + np.sqrt(-1 + 4/r) * r)
integral_proper = (np.sqrt(2)* (np.sqrt(-1 + 4/r) * r + 4 * np.arctan(np.sqrt(-1 + 4/r))))

def ef_infall():
        
    plt.semilogy(r, -1/rdot, label= r'$- 1 / \dot{r}$', c=c1, ls=':')
    plt.semilogy(r, integral_proper, label= r'Proper time [$M$]', c=c1)
    # plt.semilogy(r, vdot, label= r'$\dot{v}$')
    plt.semilogy(r, -vdot/rdot, label= r'$-\dot{v}/\dot{r}$', c=c2, ls=':')
    plt.semilogy(r, integral_coordinate, label= r'Coordinate time [$M$]', c=c2)
    
    plt.xlabel('$r / M$')
    plt.legend()

def spacetime_ef():
    plt.plot(r, integral_coordinate)

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(ef_infall)
    plot_and_save(spacetime_ef)