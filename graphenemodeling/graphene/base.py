from scipy import constants as sc
eVtoJ = sc.elementary_charge
hbar = sc.hbar

class BaseGraphene:
    """
    Base class for all types of graphene.
    Includes constants common to all types of graphene.
    Can be used to develop classes with any number of layers or any shape
    """

    def __init__(self):

        self.a = 1.42*1e-10                         # (m), Interatom spacing
        self.A = 3*(3**(1/2))*(self.a**2) / 2      # (m^2), area of unit cell of graphene
        self.g0 = 2.8*eVtoJ                         # (J), Interatom hopping potential
        self.g0prime = 0 * self.g0                 # particle hole asymmetry           
        self.vF = 3*self.a*self.g0/(2*sc.hbar)         # Fermi velocity

        self.W = 4.6 * eVtoJ                        # Work function of graphene. See 10.1038/nphys1022

        self.K = (2*sc.pi/(3*self.a)) + 1j*(2*sc.pi/(3*(3**(1/2))*self.a)) # The magnitude of the K point vector
        # Tunneling parameters
        self.kappa_elastic = 3.9*10**-10            # m^-1, inverse decay length

        # Inelastic tunneling parameters
        self.energy1 = 0.063 # V, energy of phonon
        self.sigma1 = 0.017 #  V, width of phonon peak
        self.kappa_inelastic = 2.2*10**-10 # m^-1, inverse decay length of inelastic wavefunctions
