"""
=======================================================
Constants (:mod:`graphenemodeling.graphene._constants`)
=======================================================
"""

from scipy import constants as sc
eVtoJ = sc.elementary_charge
hbar = sc.hbar

thickness =0.34e-9				# (m), thickness of graphene
a = 1.42*1e-10                         # (m), Interatom spacing
A = 3*(3**(1/2))*(a**2) / 2      # (m^2), area of unit cell of graphene
g0 = 2.8*eVtoJ                         # (J), Interatom hopping potential
g0prime = 0 * g0                 # particle hole asymmetry           
vF = 3*a*g0/(2*sc.hbar)         # Fermi velocity

W = 4.6 * eVtoJ                        # Work function of graphene. See 10.1038/nphys1022

K = (2*sc.pi/(3*a)) + 1j*(2*sc.pi/(3*(3**(1/2))*a)) # The magnitude of the K point vector

# Tunneling parameters
kappa_elastic = 3.9*10**-10            # m^-1, inverse decay length

# Inelastic tunneling parameters
energy1 = 0.063 # V, energy of phonon
sigma1 = 0.017 #  V, width of phonon peak
kappa_inelastic = 2.2*10**-10 # m^-1, inverse decay length of inelastic wavefunctions
