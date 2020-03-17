"""
=======================================================
Constants (:mod:`graphene._constants`)
=======================================================

Physical constants common in the study of graphene.
These constants are not found in ``scipy.constants``.

Carbon
======

======== =================== ====================== ==== ====
Variable Constant            Value                  Unit Ref
======== =================== ====================== ==== ====
``m_C``	 Mass of Carbon atom 1.9926468799199998e-26  kg  [1]
======== =================== ====================== ==== ====

Geometry
========
=============   ================================================================= ========= ============ ====
Variable        Constant                                                          Value     Unit         Ref
=============   ================================================================= ========= ============ ====
``thickness``   Thickness of a single layer of graphene                            3.4e-10    m          
``a``           Interatom spacing on a graphene lattice                            1.46e-10   m
``A``           Area of a unit cell of graphene lattice :math:`A=3\\sqrt{3}a^2/2`  5.24e-20   m :sup:`2`  [2]
=============   ================================================================= ========= ============ ====

Tight-Binding Parameters
========================

=========== ================================================ =========== =====
Variable    Constant                                          Value       Unit
=========== ================================================ =========== =====
``g0``		Interatom hopping potential :math:`\\gamma_0`    4.49e-19    J
``g0prime``	Particle hole asymmetry :math:`\\gamma_0'`       0           J
``vF``		Fermi velocity :math:`v_F=3\\gamma_0a/2`         9.07e5      m/s
=========== ================================================ =========== =====


Electromagnetic / Optical
=========================

=========== ==================================================== ======================= ======= ====
Variable    Constant                                             Value                   Units   Ref
=========== ==================================================== ======================= ======= ====
``W``       Work Function                                          7.37e-19                 J
``sigma_0`` Intrinsic Conductivity :math:`\\sigma_0=e^2/4\\hbar`   6.085337014469867e-05   S-m
``abs0``	Intrinsic absorption :math:`\\pi\\alpha`             0.022925309222367483    none     [3]
=========== ==================================================== ======================= ======= ====


References
==========

[1] J. Emsley, The Elements, Oxford Chemistry Guides (Oxford Univ. Press, New York, NY, 1995).

[2] Castro Neto, A.H., Guinea, F., Peres, N.M.R., Novoselov, K.S., and Geim, A.K. (2009). The electronic properties of graphene. 
Rev. Mod. Phys. 81, 109–162.
https://link.aps.org/doi/10.1103/RevModPhys.81.109.	

[3] Nair, R.R., Blake, P., Grigorenko, A.N., Novoselov, K.S., Booth, T.J., Stauber, T., Peres, N.M.R., and Geim, A.K. (2008).
Fine Structure Constant Defines Visual Transparency of Graphene. Science 320, 1308–1308.
http://www.sciencemag.org/cgi/doi/10.1126/science.1156965

"""

from scipy import constants as sc


eVtoJ = sc.elementary_charge
hbar = sc.hbar


# Carbon
m_C = 12 * sc.physical_constants['atomic mass constant'][0]

# Geometry
thickness =0.34e-9				# (m), thickness of graphene

a = 1.42*1e-10                         # (m), Interatom spacing
a1 = [a/2 * 3, a/2 * 3**0.5]
a2 = [a/2 * 3, a/2 * -(3**0.5)]
A = round(3*(3**(1/2))*(a**2) / 2 ,22)   # (m^2), area of unit cell of graphene

# Band Structure
g0 = round(2.8*eVtoJ,21)   # (J), Interatom hopping potential
g0prime = 0 * g0                 # particle hole asymmetry           
vF = round(3*a*g0/(2*sc.hbar),-3)        # Fermi velocity

W = round(4.6 * eVtoJ,21)                       # Work function of graphene. See 10.1038/nphys1022

abs0 = sc.pi*sc.alpha

K = (2*sc.pi/(3*a)) + 1j*(2*sc.pi/(3*(3**(1/2))*a)) # The magnitude of the K point vector

sigma_0	= sc.elementary_charge**2 / (4*hbar)
# Tunneling parameters
kappa_elastic = 3.9*10**-10            # m^-1, inverse decay length

# Inelastic tunneling parameters
energy1 = 0.063 # V, energy of phonon
sigma1 = 0.017 #  V, width of phonon peak
kappa_inelastic = 2.2*10**-10 # m^-1, inverse decay length of inelastic wavefunctions
