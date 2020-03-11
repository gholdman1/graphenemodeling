"""
=============================================================================
Statistical Distributions (:mod:`graphenemodeling.statistical_distributions`)
=============================================================================

These statistical distributions are not found in ``scipy.stats`` and are therefore offered here.

.. autosummary::

    FermiDirac
    BoseEinstein
    Boltzmann
    Lorentz
    Planck


"""

import numpy as np
from graphenemodeling import fundamental_constants as fc
from scipy.constants import speed_of_light,hbar, k

kB = fc.kB

def FermiDirac(E, T):
    """The Fermi-Dirac distribution.

    .. math::

        \\frac{1}{e^{E/k_B T} + 1}

    Note that we do not include a Fermi Level :math:`E_F` as this is easily implemented by exchanging :math:`E` for :math:`E-E_F`.

    Parameters
    ----------
    E:    		array-like, Energy (J)

    T:      	scalar, Temperature (K)

    Returns
    -------
    FD:         array-like, Fermi-Dirac probability of occupation of state at energy E.

    Examples
    --------

    >>> from graphenemodeling.statistical_distributions import FermiDirac
    >>> import numpy as np

    """

    # np.logaddexp reduces chance of underflow error.
    # Add a tiny offset to temperature to avoid division by zero.
    FD = np.exp( -np.logaddexp(E/(kB*(T+0.000000000001)),0) )

    return FD

def BoseEinstein(E,T):
    """
    The Bose-Einstein distribution.

    .. math::

        \\frac{1}{e^{E/k_B T} - 1}

    Note we do not include a chemical potential :math:`\\mu` as this is easily implemented by substituing :math:`E\\to E-\\mu`.

    Parameters
    ----------
    E:          array-like, Energy (J)

    T:          scalar, Temperature (K)

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

    .. math::

        e^{-E/k_B T}

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
    '''Lorentzian Response

    Not a true statistical distribution, but included here.
    
    .. math::

        A \\frac{(\\gamma/2\\pi)^2}{(x-x_0)^2 + (\\gamma/2)^2}

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
    """
    The Planck distribution.

    In terms of angular frequency :math:`\\omega`

    .. math::

        \\frac{\\hbar\\omega^3}{c^2 \\pi}\\frac{1}{e^{\\hbar\\omega/ k_B T} - 1}


    or wavelength :math:`\\lambda`

    .. math::

        \\frac{2hc^2}{\\lambda^5}\\frac{1}{e^{hc/\\lambda k_B T} - 1}

    """

    if int_var=='omega':
        
        prefactor=hbar*x**3 / (speed_of_light**2*np.pi)

        dist = 1/(np.exp(hbar*x / (k*T))-1)

    if int_var=='lambda':
        h=2*np.pi*hbar
        prefactor=2*h*speed_of_light**2 / x**5
        dist = 1/(np.exp(h*speed_of_light/(x*k*T))-1)

    return prefactor*dist
