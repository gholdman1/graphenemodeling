import numpy as np
import fundamental_constants as fc
from scipy.constants import speed_of_light,hbar, k

kB = fc.kB

def FermiDirac(E, T):
    """
    The Fermi-Dirac distribution.

    Parameters
    ----------
    E:    		array-like,
    			Energy (J)

    T:      	scalar,
    			Temperature (K)

    Returns
    ----------
    FD:         array-like,
                Fermi-Dirac probability of occupation of state at energy E.

    """

    # Using logaddexp reduces chance of underflow error
    # Adds a tiny offset to temperature to avoid division by zero.
    FD = np.exp( -np.logaddexp(E/(kB*(T+0.000000000001)),0) )

    return FD

def BoseEinstein(E,T):
    """
    The Bose-Einstein distribution.

    Parameters
    ----------
    E:          array-like,
                Energy (J)

    T:          scalar,
                Temperature (K)

    Returns
    ----------
    BE:         array-like, Bose-Einstein probability of occupation of state at energy E.

    """

    # Using logaddexp reduces chance of underflow error
    # Adds a tiny offset to temperature to avoid division by zero.
    BE = np.exp( -np.logaddexp(E/(kB*(T+0.000000000001)),1j*fc.pi) )

    return BE

def Boltzmann(E,T):
    """
    The Boltzmann distribution.

    Parameters
    ----------
    E:          array-like, Energy of state (J)

    T:          scalar, Temperature (K)

    Returns
    ----------
    boltz:      array-like, probability of occupation of state at energy E.

    """

    Tp =T + 1e-9

    #boltz = np.exp( -np.logaddexp(E/(kB*Tp),1j*fc.pi) )

    boltz = ( np.exp(E/(kB*Tp)) - 1)**(-1)

    return boltz

def Lorentz(p,x):
    '''
    Lorentzian Response

    Parameters
    ----------

    p:      length 3 array
            p[0] = response location
            p[1] = HWHM
            p[2] = response strength

    x:      array-like, points at which to evaluate Lorentzian
    '''

    return p[2] * ( (p[1]/2)/fc.pi)**2 / ( (p[0]-x)**2 + (p[1]/2)**2 )

def Planck(x,T,int_var='omega'):

    if int_var=='omega':
        
        prefactor=hbar*x**3 / (speed_of_light**2*np.pi)

        dist = 1/(np.exp(hbar*x / (k*T))-1)

    if int_var=='lambda':
        h=2*np.pi*hbar
        prefactor=2*h*speed_of_light**2 / x**5
        dist = 1/(np.exp(h*speed_of_light/(x*k*T))-1)

    return prefactor*dist