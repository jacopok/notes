from astropy.constants import codata2018 as ac
import astropy.units as u
from astropy.cosmology import Planck15 as cosmo
import numpy as np

m_u = 2.16 * u.MeV
m_d = 4.67 * u.MeV
m_n = 939.57 * u.MeV
m_p = 938.27 * u.MeV
f = 132 * u.MeV

alpha = (m_n - m_p) / (m_d - m_u) / 2
g_A = 1.26

print((m_u * m_d / (m_u + m_d) / m_n**2 * ac.hbar * ac.c).to(u.cm))

# from https://arxiv.org/pdf/hep-lat/0508009.pdf

print(((g_A * alpha / 2 / np.pi**2 / f**2 * (2 *np.log(140 * u.MeV / u.GeV)) + 1/(u.GeV)**2)* (m_u * m_d / (m_u + m_d)) * ac.hbar * ac.c).cgs)

