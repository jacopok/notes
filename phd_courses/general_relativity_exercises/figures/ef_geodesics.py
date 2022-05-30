import numpy as np
import matplotlib.pyplot as plt

CMAP = plt.get_cmap('plasma')

r = np.linspace(1e-6, 4-1e-6, num=1000)

r0 = 4
e = np.sqrt(1 - 2 / r0)

c1 = CMAP(.8)
c2 = CMAP(.2)

rdot = - np.sqrt(e**2 + 2 / r - 1)
vdot = 1 / (e - rdot)
integral_coordinate = + ((1 + np.sqrt(-1 + 4/r)) * r) + 4 * np.arctan((np.sqrt(-1 + 4/r) * (-2 + r))/(-4 + r)) + 2 * np.log(2 + np.sqrt(-1 + 4/r) * r)
integral_coordinate -= min(integral_coordinate)
integral_proper = (np.sqrt(2)* (np.sqrt(-1 + 4/r) * r + 4 * np.arctan(np.sqrt(-1 + 4/r))))

def ef_infall():
        
    plt.semilogy(r, -1/rdot, label= r'$- 1 / \dot{r}$', c=c1, ls=':')
    plt.semilogy(r, integral_proper, label= r'Proper time [$M$]', c=c1)
    # plt.semilogy(r, vdot, label= r'$\dot{v}$')
    plt.semilogy(r, -vdot/rdot, label= r'$-\dot{v}/\dot{r}$', c=c2, ls=':')
    plt.semilogy(r, integral_coordinate, label= r'Coordinate time [$M$]', c=c2)
    
    plt.xlabel('$r / M$')
    plt.legend()
    plt.ylim(1e-2, 1e2)

def spacetime_ef():

    plt.plot(r, integral_coordinate, label='Infalling massive geodesic')
    
    v_range = np.linspace(0, max(integral_coordinate), num=1000)
    
    R, V = np.meshgrid(r, v_range)
    
    dvdr_out = abs(1 / 2 / (1 - 2 / R))
    dvdr_in = np.zeros_like(dvdr_out)
    
    plt.streamplot(r, v_range, np.where(R > 2, 1, -1), dvdr_out, arrowsize=.5, density=2, linewidth=.5, color='orange')
    plt.streamplot(r, v_range, -np.ones_like(dvdr_in), dvdr_in, arrowsize=.5, density=.5, linewidth=.5, color='black')
    
    # plt.legend()
    plt.xlabel('Radius [$M$]')
    plt.ylabel('Ingoing EF coordinate $v$ [$M$]')
    plt.title('Timelike trajectory with lightlike geodesics')

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(ef_infall)
    plot_and_save(spacetime_ef)