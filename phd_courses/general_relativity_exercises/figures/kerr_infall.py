import numpy as np
import matplotlib.pyplot as plt
cmap = plt.get_cmap('plasma')

def kerr_infall():
    l = 0
    a = .9
    M = 1

    r = np.linspace(M + np.sqrt(M**2 - a**2)+1e-2, 10, num=10000)

    theta = np.pi / 2

    def metric_coeffs(r, theta=np.pi/2):
        rhosq = r**2 + (a*np.cos(theta))**2
        delta = r**2 - 2*M*r + a**2
        g_tt = - (1 - 2 * M * r / rhosq)
        g_tp = - 2 * M * a * r * np.sin(theta)**2 / rhosq
        g_pp = (r**2 + a**2 + 2 * M * r * a**2 * np.sin(theta)**2 / rhosq) * np.sin(theta)**2
        g_rr = rhosq / delta
        
        return g_tt, g_tp, g_pp, g_rr

    g_tt, g_tp, g_pp, g_rr = metric_coeffs(r[-1])
    e = np.sqrt(g_tp**2 / g_pp - g_tt)

    g_tt, g_tp, g_pp, g_rr = metric_coeffs(r)

    tdot = e / (g_tp**2 / g_pp - g_tt)
    pdot = - g_tp / g_tt * tdot

    rdotsq = ( -1 + e * tdot) / g_rr

    plt.semilogy(r, tdot, label=r'$\dot{t}$', c=cmap(.9))
    # plt.semilogy(r, abs(pdot) / tdot, label=r'$\dot{\varphi} / \dot{t}$')
    plt.semilogy(r, np.sqrt(rdotsq), label=r'$\dot{r}$', c=cmap(.5))
    
    tau = (np.cumsum(1/np.sqrt(rdotsq[:-1])[::-1]) * np.gradient(r[:-1]))[::-1]
    plt.semilogy(r[:-1], tau, label=r'$\tau$ (proper time from $r_0$)', c=cmap(.1))
    
    plt.ylabel('Various quantities, [dimensionless or in units of $M$]')
    plt.xlabel('$r$ [units of $M$]')
    plt.legend()

def kerr_trajectory():

    l = 0
    a = .1
    M = 1

    r0 = 10
    rmin = M + np.sqrt(M**2 - a**2)+1e-2
    r = rmin + np.geomspace(1e-3, r0 - rmin, num=10000)
    # r = np.linspace(rmin, r0, num=10000)
    # r = np.linspace(1, 10, num=10000)

    theta = np.pi / 2

    def metric_coeffs(r, theta=np.pi/2):
        rhosq = r**2 + (a*np.cos(theta))**2
        delta = r**2 - 2*M*r + a**2
        g_tt = - (1 - 2 * M * r / rhosq)
        g_tp = - 2 * M * a * r * np.sin(theta)**2 / rhosq
        g_pp = (r**2 + a**2 + 2 * M * r * a**2 * np.sin(theta)**2 / rhosq) * np.sin(theta)**2
        g_rr = rhosq / delta
        
        return g_tt, g_tp, g_pp, g_rr

    g_tt, g_tp, g_pp, g_rr = metric_coeffs(r[-1])
    e = np.sqrt(g_tp**2 / g_pp - g_tt)

    g_tt, g_tp, g_pp, g_rr = metric_coeffs(r)

    tdot = e / (g_tp**2 / g_pp - g_tt)
    pdot = - g_tp / g_tt * tdot

    rdotsq = ( -1 + e * tdot) / g_rr

    phi = np.cumsum(pdot / np.sqrt(rdotsq) * np.gradient(r))
    
    plt.plot(
        r * np.sin(phi),
        r * np.cos(phi),
    )
    plt.gca().set_aspect('equal')
    plt.xlabel('$x$ axis, units of $M$')
    plt.ylabel('$y$ axis, units of $M$')

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(kerr_infall)
    plot_and_save(kerr_trajectory)
