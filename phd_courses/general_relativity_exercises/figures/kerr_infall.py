import numpy as np
import matplotlib.pyplot as plt
cmap = plt.get_cmap('plasma')

M = 1

def metric_coeffs(r, theta=np.pi/2, a=.9, M=M):
    rhosq = r**2 + (a*np.cos(theta))**2
    delta = r**2 - 2*M*r + a**2
    g_tt = - (1 - 2 * M * r / rhosq)
    g_tp = - 2 * M * a * r * np.sin(theta)**2 / rhosq
    g_pp = (r**2 + a**2 + 2 * M * r * a**2 * np.sin(theta)**2 / rhosq) * np.sin(theta)**2
    g_rr = rhosq / delta
    
    return g_tt, g_tp, g_pp, g_rr

def draw_circle(radius, **kwargs):
    angles = np.linspace(0, 2*np.pi, num=1000)
    
    plt.plot(
        radius*np.cos(angles),
        radius*np.sin(angles),
        **kwargs
    )

def kerr_infall():
    l = 0
    a = .9

    r = np.linspace(0, 10, num=10000)
    # r = np.linspace(M + np.sqrt(M**2 - a**2)+1e-2, 10, num=10000)

    theta = np.pi / 2


    g_tt, g_tp, g_pp, g_rr = metric_coeffs(r[-1], a=a)
    e = np.sqrt(g_tp**2 / g_pp - g_tt)

    g_tt, g_tp, g_pp, g_rr = metric_coeffs(r, a=a)

    tdot = e / (g_tp**2 / g_pp - g_tt)
    pdot = - g_tp / g_tt * tdot

    rdotsq = ( -1 + e * tdot) / g_rr

    # plt.semilogy(r, tdot, label=r'$\dot{t}$', c=cmap(.9))
    # plt.semilogy(r, abs(pdot) / tdot, label=r'$\dot{\varphi} / \dot{t}$')
    plt.semilogy(r, np.sqrt(rdotsq), label=r'$\dot{r}$', c=cmap(.5))
    
    tau = (np.cumsum(1/np.sqrt(rdotsq[:-1])[::-1]) * np.gradient(r[:-1]))[::-1]
    plt.semilogy(r[:-1], tau, label=r'$\tau$ (proper time from $r_0$)', c=cmap(.1))
    
    plt.axvline(2*M, label='Ergosphere boundary', ls=':', c=cmap(.9))
    plt.axvline(M+np.sqrt(M**2 - a**2), label='Outer horizon', ls=':', c=cmap(.7))
    
    plt.ylabel('Various quantities, [dimensionless or in units of $M$]')
    plt.xlabel('$r$ [units of $M$]')
    plt.legend()

def kerr_trajectory():

    l = 0
    a = .35
    M = 1

    r0 = 10
    # rmin = M + np.sqrt(M**2 - a**2)+1e-2
    # r = rmin + np.geomspace(1e-3, r0 - rmin, num=10000)
    # r = np.append(
    #     2-np.geomspace(1e-3, 2-1e-3, num=1000),
    #     2+np.geomspace(1e-3, r0-2, num=1000)
    # )
    r = np.linspace(1e-4, r0, num=50_000)
    
    # r = np.linspace(rmin, r0, num=10000)
    # r = np.linspace(1, 10, num=10000)

    theta = np.pi / 2

    g_tt, g_tp, g_pp, g_rr = metric_coeffs(r[-1], a=a)
    e = np.sqrt(g_tp**2 / g_pp - g_tt)

    g_tt, g_tp, g_pp, g_rr = metric_coeffs(r, a=a)

    tdot = e / (g_tp**2 / g_pp - g_tt)
    pdot = - g_tp / g_tt * tdot

    rdotsq = ( -1 + e * tdot) / g_rr

    phi = np.cumsum(pdot / np.sqrt(rdotsq) * np.gradient(r))

    plt.scatter(
        r * np.cos(phi),
        r * np.sin(phi),
        s=.05,
        color=cmap(.1),
        label='Particle trajectory'
    )
    
    draw_circle(2*M, lw=1, c=cmap(.9), label='Ergosphere boundary')
    draw_circle(M + np.sqrt(M**2 - a**2), lw=1, c=cmap(.7), label='Outer horizon')
    
    plt.gca().set_aspect('equal')
    plt.xlabel('$x$ axis, units of $M$')
    plt.ylabel('$y$ axis, units of $M$')
    plt.legend()

def potential_barrier():

    l = 0
    a = .9

    r0 = 10
    rmin = M + np.sqrt(M**2 - a**2)+1e-2
    r = rmin + np.geomspace(1e-3, r0 - rmin, num=10000)
    # r = np.linspace(rmin, r0, num=10000)
    # r = np.linspace(1, 10, num=10000)

    theta = np.pi / 2

    g_tt, g_tp, g_pp, g_rr = metric_coeffs(r[-1])
    e = np.sqrt(g_tp**2 / g_pp - g_tt)

    g_tt, g_tp, g_pp, g_rr = metric_coeffs(r)

    tdot = e / (g_tp**2 / g_pp - g_tt)

    # plt.plot(r, e * tdot - 1)
    plt.plot(r, g_tp**2 / g_pp - g_tt)
    plt.plot(r, e**2*np.ones_like(r))

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(kerr_infall)
    plot_and_save(kerr_trajectory)
