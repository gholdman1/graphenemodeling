"""
Graphene (:mod:`graphene`)
==========================

The `graphene` module attempts to implement all types of graphene from nanostructured to planar with any number of layers.

Monolayer
---------

.. autoclass:: Monolayer
   :members:

Bilayer
-------
.. autoclass:: Bilayer
   :members:
"""


import os
import numpy as np
from scipy import integrate
from scipy import optimize
from scipy import special
from scipy import constants as sc
import warnings

from graphenemodeling import statistical_distributions as sd

from graphenemodeling import Emitter

eVtoJ = sc.elementary_charge
e0 = sc.epsilon_0
hbar = sc.hbar
Z0 = sc.physical_constants['characteristic impedance of vacuum']
sigma_0 = sc.elementary_charge**2 / (4 * sc.hbar)
kB = sc.Boltzmann




class HalfPlane(BaseGraphene):
    '''
    A class describing the properties of a half-plane of graphene.
    '''

    def __init__(self,width,type,layers):

        self.__class__= type( self.__class__.__name__,
                                (layers, BaseGraphene),
                                dict(self.__class__.__dict__))
        super(self.__class__, self).__init__()

        self.width = width

class Ribbon(BaseGraphene):
    '''
    A class describing the properties of a circular piece of graphene.
    Currently only monolayer.

    Parents
    ----------

    BaseGraphene

    layers (required argument. See Parameters below.)

    Initialization Parameters
    ----------

    width:      scalar, width of graphene ribbon.

    edge:       'armchair':
                'zig-zag':

    layers:     graphene.Monolayer for 1 layer.
                graphene.Bilayer for 2 layers

    '''

    def __init__(self,width,type,layers):

        self.__class__= type( self.__class__.__name__,
                                (layers, BaseGraphene),
                                dict(self.__class__.__dict__))
        super(self.__class__, self).__init__()

        self.width = width

class Nanostructure(BaseGraphene):
    '''
    A class describing a nanostructured piece of graphene of any shape.
    Assumed finite in all directions. Therefore, half-sheet and nanoribbons
    are not included.

    Uses the formulation as described in Ref 1 Sec. 5.3.1

    Parameters
    ---------

    eigenmodes:     n x 2 array, where n is the number of modes.
                    Pairs of normalized eigenfrequencies and mode strengths
                    Format: np.array([[freq1, stren1],[freq2,stren2],[freq3,stren3],...])

    Plan
    ----------

    Be able to pass in the parameters that describe the various modes.
    Use this class on solved systems.

    Currently only deriving properties from graphene, but could be genrealized 
    in its own file.

    Shapes:
    ----------

    Disk, Square

    Shapes that could be included
    -----------

    Pentagon, Hexagon, Triangle

    References
    ----------

    [1] Christensen 2017, Thesis
            We are using the notation of this thesis.

    [2] Miller et al. 2017, URL: http://pubs.acs.org/doi/10.1021/acs.nanolett.7b02007
    '''

    def __init__(self,layers, size=None,eigenmodes=None):

        if layers == 1:
            layerclass = Monolayer
        if layers == 2:
            layerclass = Bilayer

        # Next lines add class given in 'layers' argument
        self.__class__= type( self.__class__.__name__,
                                (layerclass, BaseGraphene   ),
                                dict(self.__class__.__dict__))
        super(self.__class__, self).__init__()

        self.size=size
        self.eigenmodes = eigenmodes

    def OpticalResponseBound(self,omega,gamma,eFermi,T):
        warnings.warn("Nanostructure.OpticalResponseBound: Only MONOLAYER is being used.")
        return Monolayer.OpticalResponseBound(self,omega,gamma,eFermi,T)

    def DipolePolarizibility(self,omega,gamma,eFermi,T,order=None):
        '''
        Equation 5.22 of Ref 1
        '''

        # Equation 5.14b of Ref [1]
        zeta = lambda w: 1j*sc.epsilon_0*1*w*(self.size) / self.OpticalConductivity(0,w,gamma,eFermi,T)

        alpha_term = lambda w, eigenmode: eigenmode[1] / (eigenmode[0] - zeta(w)) 

        alpha_sum = 0
        for em in self.eigenmodes[:order]:
            alpha_sum = alpha_sum + alpha_term(omega,em)

        alpha = 2*(self.size)**3 * alpha_sum
        
        return alpha

    def CrossSection(self,omega,gamma,eFermi,T,eps_medium,cstype):
        '''
        The cross section of a scatterer embedded in a medium.
        Equation 2.21 of Ref 1.

        Assumptions
        ----------

        Nonretarded model, medium lossless, single dipole frequency contributes.

        Parameters
        ----------

        omega:      array-like, angular frequencies of light

        eps_medium: scalar, dielectric of the medium enclosing the scatterer. Must be lossless

        cstype:     Cross Section type
                    'scattering': scattering cross section
                    'absorption': absorption cross section
                    'extinction': sum of scattering and absorption

        References
        ----------

        [1] Christensen 2017 Thesis


        '''

        alpha = self.DipolePolarizibility(omega,gamma,eFermi,T)

        kd = np.sqrt(eps_medium)*omega/sc.speed_of_light

        sigma_sca, sigma_abs = 0, 0

        if cstype == 'scattering' or cstype=='extinction':
            sigma_sca = (6*sc.pi)**(-1)*(kd**4)*np.abs(alpha)**2

        if cstype == 'absorption' or cstype=='extinction':
            sigma_abs = kd * np.imag(alpha)

        return sigma_sca + sigma_abs

class Disk(Nanostructure):
    '''
    A class describing the properties of a circular piece of graphene.
    Currently only monolayer.

    Parents
    ----------

    BaseGraphene

    layers (required argument. See Parameters below.)

    Initialization Parameters
    ----------

    radius:     scalar, radius of circle

    layers:     class, graphene.Monolayer for 1 layer.
                        graphene.Bilayer for 2 layers.

    '''

    def __init__(self,layers,radius):

        # Monolayer Graphene mode parameters
        # Normalized to L=R
        eigenfrequencies = np.array([[1.0977,2.8912],
                                     [4.9140,0.1120],
                                     [8.1337,0.0424],
                                     [11.3079,0.0224],
                                     [14.4675,0.0140],
                                     [17.6205,0.0096]])

        Nanostructure.__init__(self,layers,radius,eigenfrequencies)

    
    def DipolePolarizibility(self,omega,gamma,eFermi,T,order=None):
        return Nanostructure.DipolePolarizibility(self,omega,gamma,eFermi,T,order)

    def CrossSection(self,omega,gamma,eFermi,T,eps_medium,cstype):
        return Nanostructure.CrossSection(self,omega,gamma,eFermi,T,eps_medium,cstype)

class Ellipse(Nanostructure):
    '''
    A class describing the properties of a elliptical pieces of graphene.
    Currently only monolayer.

    Parents
    ----------

    BaseGraphene

    layers (required argument. See Parameters below.)

    Initialization Parameters
    ----------

    radius:     scalar, radius of circle

    layers:     class, graphene.Monolayer for 1 layer.
                        graphene.Bilayer for 2 layers.

    '''

    def __init__(self,layers,radius):

        # Monolayer Graphene mode parameters
        # Normalized to L=R
        eigenfrequencies = np.array([[1.0977,2.8912],
                                     [4.9140,0.1120],
                                     [8.1337,0.0424],
                                     [11.3079,0.0224],
                                     [14.4675,0.0140],
                                     [17.6205,0.0096]])

        Nanostructure.__init__(self,layers,radius,eigenfrequencies)

    
    def DipolePolarizibility(self,omega,gamma,eFermi,T,order=None):
        return Nanostructure.DipolePolarizibility(self,omega,gamma,eFermi,T,order)

    def CrossSection(self,omega,gamma,eFermi,T,eps_medium,cstype):
        return Nanostructure.CrossSection(self,omega,gamma,eFermi,T,eps_medium,cstype)
        
class Rectangle(Nanostructure):
    '''
    A class describing the properties of a circular piece of graphene.
    Currently only monolayer.

    Parents
    ----------

    BaseGraphene

    layers (required argument. See Parameters below.)

    Initialization Parameters
    ----------

    radius:     scalar, radius of circle

    layers:     class, graphene.Monolayer for 1 layer.
                        graphene.Bilayer for 2 layers.

    '''

    def __init__(self,radius,layers):
        Nanostructure.__init__(self,layers)

if __name__=="__main__":
    import doctest
    doctest.testmod(verbose=False,optionflags=doctest.ELLIPSIS)
