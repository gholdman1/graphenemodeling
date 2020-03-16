"""
=======================================================
Constants (:mod:`graphene._constants`)
=======================================================

Physical constants common to all forms of graphene.

Carbon
======

=======	===================
``m_C``	Mass of Carbon atom
=======	===================

Geometry
========

=============	=======================================
``thickness``	Thickness of a single layer of graphene
``a``			Interatom spacing on a graphene lattice
``A``			Area of a unit cell of graphene lattice
=============	=======================================

Band Structure
==============

===========	===========================
``g0``		Interatom hopping potential
``g0prime``	Particle hole asymmetry
``vF``		Fermi velocity
===========	===========================

Phonons
=======




"""

from scipy import constants as sc


eVtoJ = sc.elementary_charge
hbar = sc.hbar


# Carbon
m_C = 12 * sc.physical_constants['atomic mass constant']

# Geometry
thickness =0.34e-9				# (m), thickness of graphene
a = 1.42*1e-10                         # (m), Interatom spacing
a1 = [a/2 * 3, a/2 * 3**0.5]
a2 = [a/2 * 3, a/2 * -(3**0.5)]
A = 3*(3**(1/2))*(a**2) / 2      # (m^2), area of unit cell of graphene

# Band Structure
g0 = 2.8*eVtoJ                         # (J), Interatom hopping potential
g0prime = 0 * g0                 # particle hole asymmetry           
vF = 3*a*g0/(2*sc.hbar)         # Fermi velocity

W = 4.6 * eVtoJ                        # Work function of graphene. See 10.1038/nphys1022

K = (2*sc.pi/(3*a)) + 1j*(2*sc.pi/(3*(3**(1/2))*a)) # The magnitude of the K point vector

sigma_0	= sc.elementary_charge**2 / (4*hbar)
# Tunneling parameters
kappa_elastic = 3.9*10**-10            # m^-1, inverse decay length

# Inelastic tunneling parameters
energy1 = 0.063 # V, energy of phonon
sigma1 = 0.017 #  V, width of phonon peak
kappa_inelastic = 2.2*10**-10 # m^-1, inverse decay length of inelastic wavefunctions
