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

class BaseGraphene:
    """
    Base class for all types of graphene.
    Includes constants common to all types of graphene.
    Can be used to develop classes with more any number of layers.
    """

    def __init__(self):

        self.a = 1.42*1e-10                         # (m), Interatom spacing
        self.A = 3*np.sqrt(3)*(self.a**2) / 2      # (m^2), area of unit cell of graphene
        self.g0 = 2.8*eVtoJ                         # (J), Interatom hopping potential
        self.g0prime = 0 * self.g0                 # particle hole asymmetry           
        self.vF = 3*self.a*self.g0/(2*sc.hbar)         # Fermi velocity

        self.W = 4.6 * eVtoJ                        # Work function of graphene. See 10.1038/nphys1022

        self.K = (2*sc.pi/(3*self.a)) + 1j*(2*sc.pi/(3*np.sqrt(3)*self.a)) # The magnitude of the K point vector
        # Tunneling parameters
        self.kappa_elastic = 3.9*10**-10            # m^-1, inverse decay length

        # Inelastic tunneling parameters
        self.energy1 = 0.063 # V, energy of phonon
        self.sigma1 = 0.017 #  V, width of phonon peak
        self.kappa_inelastic = 2.2*10**-10 # m^-1, inverse decay length of inelastic wavefunctions


class Monolayer(BaseGraphene):

    def __init__(self,mobility=None,thickness=0.34e-9):
        '''
        Parameters
        ----------

        mobility:   tuple, [mobility, Tmeas] of mobility (m^2/V-s)
                            and temperature at which it was taken
        '''
        BaseGraphene.__init__(self)

        if mobility != None:
            self.mu0 = mobility[0]
            self.mu0T= mobility[1]

        self.thickness = thickness

        
    ##########
    # Basics #
    ##########

    def Hamiltonian(self,k):
        '''
        Tight-binding Hamiltonian in k-space.

        Parameters
        ----------

        k:      array-like, wavevector of carrier. Use complex k=kx + iky for 2D wavevectors.

        Returns
        ----------

        H:      2x2 array, Tight-binding Hamiltonian evaluated at k.

        References
        ----------

        [1]     Castro Neto et al. Reviews of Modern Physics 81, 2009.
                URL:

        '''

        pass

    def DiracFermionDispersion(self,k,model,eh=1):
        '''
        Gives the dispersion of Dirac fermions in monolayer graphene.

        Parameters
        ----------

        k:          array-like, wavevector of Dirac fermion relative to K vector
                                For 2D wavevectors, use comlex k= kx + i ky


        model:      'LowEnergy': Linear approximation of dispersion.
                    'FullTightBinding': Eigenvalues of tight-binding approximation. 

        eh:         1 or -1, return electrons or holes respectively

        Returns
        ----------

        energy:     array-like, energy of Dirac fermion with wavenumber k.

        References
        ----------

        [1] 

        Examples
        --------
        Plot the Fermion dispersion relation.

        >>> from graphenemodeling import graphene
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> mlg = graphene.Monolayer()
        >>> from scipy.constants import elementary_charge as eV
        >>> eF = 0.4*eV
        >>> kF = mlg.kFermi(eF,model='LowEnergy')
        >>> k = np.linspace(-2*kF,2*kF,num=100)
        >>> conduction_band = mlg.DiracFermionDispersion(k,model='LowEnergy')
        >>> valence_band = -conduction_band
        >>> fig, ax = plt.subplots(figsize=(5,6))
        >>> ax.plot(k/kF,conduction_band/eF,'k')
        [...
        >>> ax.plot(k/kF,valence_band/eF, 'k')
        [...
        >>> ax.plot(k/kF,np.zeros_like(k),color='gray')
        [...
        >>> ax.axvline(x=0,ymin=0,ymax=1,color='gray')
        <...
        >>> ax.set_axis_off()
        >>> plt.show()

        Plot the full multi-dimensional dispersion relation.

        >>> from graphenemodeling import graphene
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> from mpl_toolkits import mplot3d # 3D plotting
        >>> mlg = graphene.Monolayer()
        >>> kmax = np.abs(mlg.K)
        >>> emax = mlg.DiracFermionDispersion(0,model='FullTightBinding')
        >>> kx = np.linspace(-kmax,kmax,num=100)
        >>> ky = np.copy(kx)
        >>> k = (kx + 1j*ky[:,np.newaxis]) + mlg.K # k is relative to K. Add K to move to center of Brillouin zone
        >>> conduction_band = mlg.DiracFermionDispersion(k,model='FullTightBinding',eh=1)
        >>> valence_band = mlg.DiracFermionDispersion(k,model='FullTightBinding',eh=-1)
        >>> fig = plt.figure(figsize=(8,8))
        >>> fullax = plt.axes(projection='3d')
        >>> fullax.view_init(20,35)
        >>> KX, KY = np.meshgrid(kx,ky)
        >>> fullax.plot_surface(KX/kmax,KY/kmax,conduction_band/mlg.g0,rstride=1,cstride=1,cmap='viridis',edgecolor='none')
        <...
        >>> fullax.plot_surface(KX/kmax,KY/kmax,valence_band/mlg.g0,rstride=1,cstride=1,cmap='viridis',edgecolor='none')
        <...
        >>> fullax.set_xlabel('$k_x/|K|$')
        Text...
        >>> fullax.set_ylabel('$k_y/|K|$')
        Text...
        >>> fullax.set_zlabel('$\\epsilon/\\gamma_0$')
        Text...
        >>> fullax.set_title('Brillouin Zone of Graphene')
        Text...
        >>> plt.show()
        '''


        if model == 'LowEnergy':

            return eh*sc.hbar*self.vF*np.abs(k)

        if model == 'FullTightBinding':

            k = k - self.K
            f = lambda k: (2*np.cos(np.sqrt(3)*np.imag(k)*self.a) 
                            + 4*np.cos((np.sqrt(3)*np.imag(k)/2)*self.a)*np.cos((3/2)*np.real(k)*self.a) )

            return eh*self.g0*np.sqrt(3+ f(k)) - self.g0prime*f(k)

    def kFermi(self,eFermi,model):
        '''
        The Fermi wavevector, i.e. the wavevector of the state at
        the Fermi energy.

        Parameters
        ----------

        eFermi:      array-like, Fermi level

        model:       'LowEnergy': Value gets derived from linear approximation of dispersion.
 
        Examples
        --------
        Confirm energy of Fermi wavevector is equal to Fermi level.

        >>> from graphenemodeling import graphene
        >>> from scipy.constants import elementary_charge as eV
        >>> mlg = graphene.Monolayer()
        >>> eF = 0.4 * eV # Fermi level is 0.4 eV
        >>> kF = mlg.kFermi(eF, model='LowEnergy')
        >>> mlg.DiracFermionDispersion(kF,model='LowEnergy')/eV
        0.4
        '''

        if model == 'LowEnergy':
            return np.abs(eFermi) / (sc.hbar*self.vF)

        if model == 'FullTightBinding':
            '''
            Code to numerically invert DiracFermionDispersion(kF,model='FullTightBinding')

            Likely would use a root-finding procedure
            '''
            pass

    def DensityOfStates(self,E,model):
        '''
        The density of states per unit cell in graphene at energy E.

        Parameters
        ----------

        model:      'LowEnergy': DOS derived from linear approximation of dispersion.

        References
        ----------

        [1]     Castro Neto et al. Reviews of Modern Physics 81, 2009.
                URL:

        Examples
        --------
        Plot the density of states in the `LowEnergy` approximation and `FullTightBinding` model.

        >>> from graphenemodeling import graphene
        >>> import matplotlib.pyplot as plt
        >>> mlg = graphene.Monolayer()
        >>> E = np.linspace(-3,3,num=200) * mlg.g0
        >>> DOS_low = mlg.DensityOfStates(E,model='LowEnergy')
        >>> DOS_full = mlg.DensityOfStates(E,model='FullTightBinding')
        >>> plt.plot(E/mlg.g0,DOS_full/np.max(DOS_full),label='FullTightBinding')
        [<...
        >>> plt.plot(E/mlg.g0,DOS_low/np.max(DOS_full),label='LowEnergy')
        [<...
        >>> plt.xlabel('$E/g_0$')
        Text...
        >>> plt.ylabel('DOS (a.u.)')
        Text...
        >>> plt.legend()
        <...
        >>> plt.show()
        '''

        if model=='LowEnergy':

            E = np.abs(E)

            return 2 * E / (sc.pi*(sc.hbar*self.vF)**2)

        elif model=='FullTightBinding':
            '''
            Equation 14 of Ref 1
            '''

            prefactor = 4*np.abs(E) / (sc.pi*self.g0)**2

            def fZ0(E):
                if np.abs(E)<np.abs(self.g0):
                    term1 = (1+np.abs(E/self.g0))**2
                    term2 = -((E/self.g0)**2 - 1)**2 / 4

                    return term1 + term2

                else: return 4*np.abs(E/self.g0)

            def fZ1(E):
                if np.abs(E)<np.abs(self.g0):
                    return 4*np.abs(E/self.g0)

                else:
                    term1 = (1+np.abs(E/self.g0))**2
                    term2 = -((E/self.g0)**2 - 1)**2 / 4

                    return term1 + term2

            dos = np.empty_like(E)

            for j, e in np.ndenumerate(E):
                Z0 = fZ0(e)
                Z1 = fZ1(e)
                ellip = special.ellipk(np.sqrt(Z1/Z0))

                dos[j] = (1/np.sqrt(Z0)) * ellip
            return prefactor * dos /self.A
        else:
            print('The model %s is not available' % (model))

    def CarrierDensity(self,muC,T,model):
        '''
        Computes the carrier density at nonzero temperature directly from the band structure.
        '''

        warnings.warn('Monolayer.CarrierDensity has not been verified')
        n = np.empty_like(muC)

        for i, m in np.ndenumerate(muC):
            p_electron = lambda e: self.DensityOfStates(e,model) * sd.FermiDirac(e-m,T)
            p_hole = lambda e: self.DensityOfStates(e,model) * (1 - sd.FermiDirac(e-m,T))
            n[i] = ( integrate.quad(p_electron,0,3*self.g0,points=(self.g0,m))[0] -
                     integrate.quad(p_hole,-3*self.g0,0,points=(-self.g0,-m))[0]   )
        return n

    def FermiLevel(self,n,T=0):
        '''
        Returns the Fermi level given the carrier density.
        '''

        if T==0:
            return sc.hbar*self.vF*np.sqrt(sc.pi*n)

        else:
            warnings.warn('Monolayer.FermiLevel: not available')

    ########################
    # Electrical Transport #
    ########################
    
    def Mobility(self,Te,mu0,T0):
        '''
        Temperature dependent mobility.

        See page 4 of Reference 1.

        35, 37, 38 

        References
        ----------

        [1] Shiue et al. 2019
            URL: http://www.nature.com/articles/s41467-018-08047-3

        [2] Dorgan et al. 2013.
            URL: https://doi.org/10.1021/nl400197w

        [3] Bae et al. 2010
            URL: https://doi.org/10.1021/nl1011596
        '''

        mu = mu0*(Te/T0)**2.3

        return mu

    def ScatteringRate(self,mobility,eFermi):
        '''
        Estimated DC scattering rate from mobility.

        Paremeters
        ----------

        mobility:   scalar, mobility (m^2/V-s)

        eFermi:     scalar, Fermi level (J)

        Returns
        ----------

        rate:       scalar, scattering rate
        '''

        # Scattering time
        tau = mobility*eFermi / (sc.elementary_charge*self.vF**2)

        rate = 1/tau
        return rate

    def Polarizibility(self,q,omega,gamma,eFermi,T=0):
        '''
        The Polarizibiliy function in graphene. Can be derived from a
        self-consistent field method or the Kubo formula.

        For gamma == 0, this returns equation 17 of Ref 2, which is the
        polarizibility for general complex frequencies.

        For gamma != 0, we return the Mermin-corrected Relaxation time approximation
        (Eqn 4.9 of Ref 1), which calls the gamma==0 case.

        Parameters
        ----------

        q:      array-like, difference between scattered wavevector and incident

        omega:  array-like, frequency

        gamma:  scalar, scattering rate due to mechanisms such as impurities (i.e. not Landau Damping)
                        We use the Mermin-corrected Relaxation time approximation (Eqn 4.9 of Ref 1)

        eFermi: scalar, Fermi level of graphene.

        References
        ----------

        [1] Christensen Thesis 2017

        [2] Sernelius 2012

        [2] Wunsch 2006

        [3] Hwang and Das Sarma 2007

        '''

        if gamma==0 and T==0:
            '''
            Equation 17 of Ref 2
            '''

            prefactor = -self.DensityOfStates(eFermi,model='LowEnergy')

            kF = self.kFermi(eFermi, model='LowEnergy')

            x = q / (2*kF)
            zbar = sc.hbar*omega / (2*eFermi)

            f = lambda x,zbar: (np.arcsin( (1-zbar)/x) + np.arcsin( (1+zbar)/x )
                            - ((zbar-1)/x)*np.sqrt(1 - ((zbar-1)/x)**2 )
                            + ((zbar+1)/x)*np.sqrt(1 - ((zbar+1)/x)**2 ) )


            dd = 1 + (x**2 / (4*np.sqrt(x**2 - (zbar+1e-9*1j)**2 ))) * (sc.pi - f(x,zbar+1e-9*1j))

            return prefactor * dd

        elif gamma !=0:
            # Mermin-corrected Relaxation Time Approximation (Eqn 4.9 of Ref 1)
            pol_complex_arg = self.Polarizibility(q,omega+1j*gamma,0,eFermi,T=0)
            pol_0 = self.Polarizibility(q,0,0,eFermi,T=0)

            numerator = (1 + 1j*gamma/omega) * pol_complex_arg
            denominator = 1 + ( 1j*gamma/omega * pol_complex_arg / pol_0 )
            return numerator / denominator

    def dPolarizibility(self,q,omega,gamma,eFermi,T,dvar,diff=1e-7):
        '''
        Returns the derivative of the real part of the polarizibility at q, omega
        with respect to the chosen variable, dvar.

        Parameters
        ----------

        q:      array-like,

        omega:

        gamma:  scalar, the scattering rate in units (1/s)

        eFermi: scalar, the Fermi level (J)

        T:      scalar, Temperature (K)

        dvar:   'omega': Take the partial wrt omega
                'q': Take the partial wrt q

        diff:   Size of finite different to use when computing the derivative.
                Method uses central difference.
        '''

        if dvar == 'omega':
            P = lambda w: np.real(self.Polarizibility(q,w,gamma,eFermi,T))
            wa, wb = omega*(1-diff), omega*(1+diff)
            return (P(wb)-P(wa))/(2*omega*diff)

        elif dvar == 'q':
            P = lambda qv: np.real(self.Polarizibility(qv,omega,gamma,eFermi,T))
            qa,qb = q*(1-diff), q*(1+diff)
            return (P(qb)-P(qa))/(2*q*diff)

    ######################
    # Optical Properties #
    ######################

    def OpticalConductivity(self,q,omega,gamma,eFermi,T,model=None):
        '''
        Return the diagonal conductivity of graphene sigma_xx.

        Parameters
        ----------

        q:          array-like, wavenumbers at which to evaluate the nonlocal conductivity.
                                Choose q=0 for the LOCAL conductivity.

        omega:      array-like, frequency

        eFermi:     scalar, the Fermi energy (J)

        gamma:      scalar, scattering rate

        T:          scalar, Temperature

        model:      Typically 'None', but for a specific model, specify it.


        Returns
        ----------

        conductivity:   array-like, conductivity at every value of omega

        References:
        ----------


        '''

        # Local case
        if np.all(q) == 0:

            if T!=0:
                intra_pre = 4 * sigma_0 * (2*1j*kB*T) / (sc.pi*sc.hbar)
                inter_pre = sigma_0

                ### Intraband Contribution ###

                # Using np.logaddexp() avoids the large cosh in ln( cosh(1/T) )
                x = eFermi / (2*kB*T)
                intra = lambda w: intra_pre * ( 1 / (w + 1j*gamma) ) * np.logaddexp(x,-x)

                ### Interband Contribution ###

                H = lambda energy: sd.FermiDirac(-energy-eFermi,T) - sd.FermiDirac(energy-eFermi,T)

                integrand = lambda energy,w: ( H(energy) - H(sc.hbar*w/2) ) / (sc.hbar**2 * (w + 1j*gamma)**2 - 4 * energy**2)

                def integral(w):
                    result = np.empty_like(w,dtype=complex)

                    for i, frequency in np.ndenumerate(w):
                        integrand_re = lambda e: np.real(integrand(e,frequency))
                        integrand_im = lambda e: np.imag(integrand(e,frequency))

                        result[i] =(     integrate.quad(integrand_re,0,10*eFermi,points=(eFermi/sc.hbar,2*eFermi/sc.hbar))[0]
                                    + 1j*integrate.quad(integrand_im,0,10*eFermi,points=(eFermi/sc.hbar,2*eFermi/sc.hbar))[0] )

                    return result

                inter = lambda w: inter_pre * ( H(sc.hbar * w / 2) +
                                                (4*1j/sc.pi) * sc.hbar*(w + 1j*gamma)*integral(w) )

                conductivity= intra(omega) + inter(omega)

            if T==0:
                intra = lambda w: 1j*sigma_0 * 4*eFermi / (sc.pi*sc.hbar* (w + 1j*gamma))

                inter = lambda w: sigma_0 * ( np.heaviside(sc.hbar*w - 2*eFermi,0.5) + 
                                                (1j/sc.pi) * np.log(np.abs((2*eFermi-sc.hbar*w)/(2*eFermi+sc.hbar*w))))

                conductivity = intra(omega) + inter(omega)

        # Nonlocal case
        else:
            if T==0:
                conductivity = 1j*sc.elementary_charge**2 * (omega / q**2) * self.Polarizibility(q,omega,gamma,eFermi,T)

            else:
                pass

        return conductivity

    def OpticalConductivityMatrix(self,q,omega,gamma, eFermi,T):
        '''
        Returns the conductivity matrix of monolayer graphene.

        Parameters
        ----------

        omega:

        Returns
        ----------

        sigma_matrix:   2x2 numpy array, conductivity of monolayer graphene

        '''


        conductivity_matrix = np.array([[self.OpticalConductivity(q,omega,gamma,eFermi,T),0],
                                        [0,self.OpticalConductivity(q,omega,gamma,eFermi,T)]])

        return conductivity_matrix

    def Permittivity(self,omega,eFermi,T, gamma=None,epsR=None,model=None):
        '''
        Returns the in-plae permittivity of graphene.

        Assumes local conductivity

        Parameters
        ----------

        epsR:       scalar, background relative permittivity

        References
        ----------

        [1] “Lumerical: Modeling Methodology.” n.d. Accessed April 1, 2019.
            https://apps.lumerical.com/other_application_graphene_simulation_tips.html.

        '''

        if model=='Lumerical:Falkovsky':
            '''
            Use eqn 10 of Ref 1
            '''
            x1 = sc.elementary_charge
            x2 = sc.hbar
            x3 = eFermi
            x4 = self.vF
            x5 = self.Mobility(T,self.mu0,self.mu0T) # mobility at the new temperature
            x6 = epsR

            x7 = self.thickness # 3.4 angstroms by default
            x8 = sc.epsilon_0

            prefactor = x1**2*x3 / ( sc.pi * x2**2)
            denominator = omega**2 + ( x1*x4**2 / (x5*x3) )**2

            term1 = - prefactor*(omega*x8*x7)**(-1) * omega / denominator
            term2 = 1j*prefactor * (x1*x4**2 / (omega*x5*x3*x8*x7)) / denominator

            eps = x6 + term1 + term2

            return eps

        else:
            eps = 1 + 1j*self.OpticalConductivity(0,omega,gamma,eFermi,T)/(omega*sc.epsilon_0)

        return eps

    def FresnelReflection(self,kpar,omega,gamma,eFermi,T,eps1,eps2,polarization):
        '''
        The Fresnel Reflection coefficients of light incident from above (medium 1, eps1).

        Equation 5.4 of Ref 1
        
        Parameters
        ----------

        kp:             Parallel (in-plane) momentum of indicent light.

        omega:          Frequency of incident light.

        eps1, eps2:     Permittivities above and below graphene, respectively.
                        Could also be made callable such that eps1=eps1(kpar,omega)

        polarization:   's'/'TE' or 'p'/'TM' for s- or p-polarization.

        References
        ----------

        [1] Christensen Thesis 2017
        '''

        kperp1, kperp2 = np.sqrt(eps1*(omega/sc.speed_of_light)**2 - kpar**2 + 1e-9*1j), np.sqrt(eps2*(omega/sc.speed_of_light)**2 - kpar**2 + 1e-9*1j)

        sigma = self.OpticalConductivity(kpar,omega,gamma,eFermi,T)

        if polarization=='p' or polarization=='TM':
            numerator   = eps2*kperp1 - eps1*kperp2 + ( sigma*kperp1*kperp2 / (sc.epsilon_0*omega) )
            denominator = eps2*kperp1 + eps1*kperp2 + ( sigma*kperp1*kperp2 / (sc.epsilon_0*omega) )

        if polarization=='s' or polarization=='TE':
            numerator = kperp1 - kperp2 - sc.mu_0*omega*sigma
            denominator = kperp1 + kperp2 + sc.mu_0*omega*sigma

        return numerator / denominator

    def LocalDensityOfStates(self,d,omega,gamma,eFermi,T):
        '''
        The LDOS a distance d above a plane of graphene.

        Eqn 44 of SM of Ref 1

        Assuning plane in vauum for now.

        References
        -----------

        [1] Miller et al. 2017

        '''

        k0 = (omega/sc.speed_of_light)
        ldos0 = k0**2 / (2*sc.pi**2*sc.speed_of_light) # Free space LDOS

        integral = np.empty_like(d)

        for i, d0 in np.ndenumerate(d):
            k0w = (omega/sc.speed_of_light)
            Im_rp      = lambda kpar: np.imag( self.FresnelReflection(kpar,omega,gamma,eFermi,T,1,1,'p') )

            integrand = lambda kpar: (kpar**2/k0w**3)*Im_rp(kpar)*np.exp(-2*kpar*d0)

            a,b = 1e-3, np.abs(self.K) # Increasing this bound does not lead to more convergence

            k_plasmon=self.InversePlasmonDispersion(omega,gamma,eFermi,1,1,T,model='nonlocal')

            integral[i] = integrate.quad(integrand,a,b,
                                        points=(k_plasmon),limit=100)[0]

        return ldos0 * integral      

    def OpticalResponseBound(self,omega,gamma,eFermi,T,d=None,method='svd',restype='material'):
        '''
        A general material optical response bound independent of geometrical parameters.

        See equation 6 of Ref 1 (and remove the beta).

        Assumes a linear conductivity model K=sigma * E, could be generalized
        to K = L E where L is a differential operator.

        For LDOS options, assumes a planar geometry. See Ref 1 for solid angle correction for nonplanar

        Parameters
        ----------

        omega:      array-like, the frequency range of interest (rad/s)

        gamma:      scalar, the scattering rate (rad/s)

        eFermi:     scalar, Fermi level (J)

        T:          scalar, Temperature (K)

        d:          scalar, Distance from graphene (m), if using LDOS restype (see below)

        method:     'svd': Use singular value decomposition version of formula
                    'scalar': Uses the simplified scalar version.

        restype:    Bound on response types given below
                    'material': Instrinsic material FOM
                    'CSabs': Absorption cross section / area
                    'CSext': Extinction cross section / area
                    'CSsca': Scattering cross section / area
                    'LDOStot': Total LDOS / Free space LDOS
                    'LDOSnrad': Nonradiative LDOS / Free space LDOS
                    'LDOSrad': Radiative LDOS / Free space LDOS


        References
        ----------

        [1] Miller et al. 2017
            URL: http://pubs.acs.org/doi/10.1021/acs.nanolett.7b02007
        '''

        if np.all(d):
            k0 = omega/sc.speed_of_light
            LDOSprop = 3 / (8 * (k0*d)**4)

        prop = dict({'material':1,
                    'CSabs':1, 'CSext': 1, 'CSsca': 1/4,
                    'LDOStot':1*LDOSprop,'LDOSnrad':1*LDOSprop,'LDOSrad':LDOSprop/4})

        if method == 'scalar':
            sigma = self.OpticalConductivity(0,omega,gamma,eFermi,T)

            bound = Z_0 * np.abs(sigma)**2 / np.real(sigma)

        elif method == 'svd':
            bound = np.empty_like(omega)

            for i, w in np.ndenumerate(omega):
                s = self.OpticalConductivityMatrix(0,w,gamma,eFermi,T)
                s_dag = np.conj(np.transpose(s))
                s_re_inv = np.linalg.inv(np.real(s))
                prod = np.dot(s_dag,np.dot(s_re_inv,s))
                two_norm= np.linalg.svd(prod)[1][0]
                bound[i] = Z_0*two_norm



        return prop[restype] * bound

    ###########
    # Phonons #
    ###########

    def PhononSelfEnergy(self):
        pass

    ##############
    # Plasmonics #
    ##############

    def PlasmonDispersion(self,q,gamma,eFermi,eps1,eps2,T,model):
        '''
        Returns the nonretarded plasmon dispersion relations E(q) for a surface
        plasmon in an infinite sheet of graphene sandwiched between two
        dielectrics eps1 and eps2.

        All values returned are assumed to be at zero temperature with no loss (gamma).

        Parameters
        ----------

        q:          array-like, wavenumber of the plasmon

        eps1,eps2:  scalar, the relative permittivity of each dielectric

        model:      'intra':    intraband dispersion
                    'local':    uses the intraband + interband constributions
                                to the conductivity.
                    'nonlocal': Uses fully nonlocal conductivity to get dispersion

        Returns
        ----------

        omega:      array-like, the frequency of the plasmon with wavenumber q


        '''

        epsavg = (eps1+eps2)/2

        # Analytical expression in intraband approximation
        if model=='intra':
            radical = q * sc.elementary_charge**2 * eFermi / (2*sc.pi * sc.epsilon_0 * epsavg)
            return (1/sc.hbar) * np.sqrt(radical)

        if model=='local':
            omega = np.empty_like(q)

            sigma = lambda w: self.OpticalConductivity(0,w,gamma,eFermi,T=0)

            for i,q0 in np.ndenumerate(q):
                root_eqn = lambda w: 1 - np.imag(sigma(w))*q0 / (2*sc.epsilon_0*epsavg*w)

                a, b = self.PlasmonDispersion(q0,gamma,eFermi,eps1,eps2,T,model='intra'), 1e-8
                omega[i] = optimize.brentq(root_eqn,a,b)

            return omega

        if model=='nonlocal':
            omega = np.empty_like(q)

            kF = self.kFermi(eFermi,model='LowEnergy')

            for i, q0 in np.ndenumerate(q):
                root_eqn = lambda w: self.PlasmonDispersionRoot(q0,w,gamma,eFermi,eps1,eps2,T=0)
                
                # Frequency is definitely below 1,1 intraband dispersion
                b = self.PlasmonDispersion(q0,gamma,eFermi,1,1,T,model='intra')

                # Frequency is definitely above the minimum which should be <0
                a = optimize.minimize_scalar(root_eqn,bounds=((eFermi/sc.hbar)*q0/kF,b),method='bounded').x
                
                if root_eqn(a) > 0:
                    #warnings.warn("Monolayer.PlasmonDispersion(model='nonlocal'): No root exists for q=%.5f" % (q0))
                    omega[i]=0
                else:
                    root_eqn_abs= lambda w: np.abs(root_eqn(w))

                    omega[i] = optimize.minimize_scalar(root_eqn_abs,bounds=(a,b),method='bounded').x

            # Maybe add in a fit to the width for the nonlocal version

            return omega

    def PlasmonDispersionRes(self,q,gamma,eFermi,eps1,eps2,T,exp_res=1):
        '''
        Uses the FresnelReflection coefficients to numerically search for the 
        plasmon dispersion

        Parameters
        ----------

        exp_res:    expected number of resonances
        '''
        q = np.atleast_1d(q)
        ErrorFunc = lambda p,x,y: sd.Lorentz(p,x) - y

        pFit = np.empty((np.size(q),3))

        for i,q0 in np.ndenumerate(q):
            w1 = q0*self.vF
            w2 = self.PlasmonDispersion(q0,gamma,eFermi,eps1,eps2,T,model='intra')
            w0 = (w2+w1)/2
            # omega=np.linspace(q0*self.vF,2*q0*self.vF,num=300)
            omega=np.linspace(w1,w2,num=300)

            y = np.imag(self.FresnelReflection(q0,omega,gamma,eFermi,T,eps1,eps2,'TM'))

            p0=[w0,0.1*w0,10]

            fit = optimize.leastsq(ErrorFunc,p0,
                                        args=(omega,y))

            pFit[i,:] = fit[0]



        return np.abs(pFit)

    def InversePlasmonDispersion(self,omega,gamma,eFermi,eps1,eps2,T,model):
        '''
        Returns the wavenumber of a plasmon given the frequency.

        Useful when determining the wavelength of a plasmon excited by light.
        '''

        kF = self.kFermi(eFermi,model='LowEnergy')

        cutoff = 4*kF
        q = np.empty_like(omega)

        for i, omega in np.ndenumerate(omega):
            root_eqn = lambda q: np.abs( omega - self.PlasmonDispersion(q,gamma,eFermi,eps1,eps2,T,model) )

            reps = 1
            while reps < 5:
                q[i] = optimize.minimize_scalar(root_eqn,bounds=(1e-6,reps*cutoff),method='bounded').x

                if q[i] >= cutoff:
                    reps=reps+1
                else:
                    reps=5

        return q

    def PlasmonDispersionRoot(self,q,omega,gamma,eFermi, eps1,eps2 ,T):
        '''
        The equation used for numerically solving the plasmon dispersion in the nonretarded regime.
        '''

        epsavg = (eps1+eps2)/2

        return 1 - np.imag(self.OpticalConductivity(q,omega,gamma,eFermi,T))*q / (2*sc.epsilon_0*epsavg*omega)

    def PlasmonDispersionRelation(self,q,omega,gamma,eFermi,eps1,eps2,T):
        '''
        Computes the Full Dispersion relation.

        Parameters
        ----------

        q:      array-like (complex), 

        omega:  array-like (complex),
        '''

        
        pass

    def PlasmonDispersionLoss(self,omega,gamma,eFermi,eps1,eps2,T,model):
        '''
        The loss of the plasmon wavenumber q=q1+iq2. Returns q2. Equation 15 of Ref [1]
        with tau = infinity (or gamma = 0). Assumes q2<<q1

        References
        ----------

        [1] Jablan et al. 2009 (note tau = 1/ gamma)
        '''
        
        if model == 'nonlocal':
            q1 = self.InversePlasmonDispersion(omega,gamma,eFermi,eps1,eps2,T,model)

            pol = self.Polarizibility(q1,omega,gamma,eFermi,T=0)
            pol0= self.Polarizibility(q1,1e-9,gamma,eFermi,T=0)

            dpolq = self.dPolarizibility(q1,omega,gamma,eFermi,T,dvar='q')
            dpolw = self.dPolarizibility(q1,omega,gamma,eFermi,T,dvar='omega')

            numerator = np.imag(-pol) + gamma * (-1)*dpolw + (gamma/omega) * np.real(-pol * (1- (pol/pol0)))
            denominator = (1/q1) * np.real(-pol) - (-1)*dpolq

            q2 = numerator / denominator

        return q2

    def dPlasmonDispersion(self,q,gamma,eFermi,eps1,eps2,T,model,dvar=None,diff=1e-7):
        '''
        Derivative of plasmon frequency with respect to dvar.

        Parameters
        ----------

        q:      array-like,

        omega:

        gamma:  scalar, the scattering rate in units (1/s)

        eFermi: scalar, the Fermi level (J)

        T:      scalar, Temperature (K)

        dvar:   'omega': Take the partial wrt omega
                'q': Take the partial wrt q

        diff:   Size of finite different to use when computing the derivative.
                Method uses central difference.
        '''

        if dvar == 'eFermi':
            result = np.empty_like(q)
            for i, q0 in np.ndenumerate(q):
                w = lambda eF: self.PlasmonDispersion(q0,gamma,eF,eps1,eps2,T,model)
                e1, e2 = eFermi*(1-diff), eFermi*(1+diff)
                result[i] = (w(e2)-w(e1))/(2*eFermi*diff)

            return result

        if dvar=='q':
            pass

        pass

    # def DipoleDecayRate(self,z,omega,gamma,eFermi,T,eps1,eps2):
    #     '''
    #     Decay rate of a dipole emitter placed near graphene.
    #     Right now, the dipole only points perpendicular.

    #     This should be moved to the Emitter.Dipole object

    #     Eqn 5 in the SM of Ref 1.

    #     References
    #     ----------

    #     [1] Koppens et al. 2011
    #         URL: https://doi.org/10.1021/nl201771h
    #     '''
    #     warnings.warn('Monolayer.DipoleDecayRate: Ignoring dipole components in xy-plane')

    #     dipole = Emitter.Dipole()

    #     d = np.array([0,0,1])

    #     integral = np.empty_like(omega)

    #     for i, w in np.ndenumerate(omega):

    #         kperp   = lambda kpar: np.sqrt(eps1*(w/sc.speed_of_light)**2 - kpar**2)
    #         rp      = lambda kpar: self.FresnelReflection(kpar,w,gamma,eFermi,T,eps1,eps2,'p')
    #         perpterm = lambda kpar: 2*np.abs(d[2])**2 * kpar**2 * rp(kpar)

    #         integrand = lambda kpar: np.real( (kpar - 1e-9*1j) * np.real( perpterm(kpar-1e-9*1j) * np.exp(2*1j*kperp(kpar-1e-9*1j)*z) / kperp(kpar-1e-9*1j) ) )

    #         b = np.abs(self.K) # Increasing this bound does not lead to more convergence
    #         kpar_pol = np.sqrt(eps1) * (w/sc.speed_of_light)
    #         k_plasmon=self.InversePlasmonDispersion(w,gamma,eFermi,eps1,eps2,T,model='nonlocal')

    #         integral[i] = integrate.quad(integrand,1e-3,b,
    #                                     points=(kpar_pol,k_plasmon),limit=100)[0]

    #     return dipole.DecayRate(omega,d) + (1/sc.hbar) * integral


class Bilayer(BaseGraphene):
    """
    Bilayer graphene class which inherits the parameters of the
    BaseGraphene class.
    """

    def __init__(self):
        BaseGraphene.__init__(self)
        g1  = 0.358 * eVtoJ # (J), A1-B1 hopping potential
        g3  = 0.3   * eVtoJ # (J), A1-B2 hopping potential
        g4  = 0.12  * eVtoJ # (J), A1-A2 hopping potential (McCann Koshino 2013)
        d   = 3*(10**-10)   # (m), interlayer spacing
        approx_choices = ['None', 'Common', 'LowEnergy']
        C = e0 / d

    def Hamiltonian(self,k,u):
        '''
        Returns the full tight-binding Hamiltonian of BLG.
        For array-like inputs of k, the Hamiltonian of the
        ith value of k is Hamiltonian[:,:,i]

        Parameters
        ----------
        k:  array-like
            Wavenumber (1/m).

        u:  scalar
            Interlayer potential energy difference (J).

        Returns
        ----------
        H:  array-like
            Tight-binding Hamiltonian of bilayer graphene.
        '''

        k = np.atleast_1d(k)
        length = np.shape(k)[0]
        ones = np.ones(length)

        # Diagonals
        H11 = H22 = -u/2 * ones
        H33 = H44 =  u/2 * ones

        # Intralayer
        H12 = H21 = H34 = H43 = sc.hbar * self.vF * k

        # Interlayer A1-B2
        H23 = H32 = self.g1 * ones

        # Trigonal Warping
        H14 = H41 = np.sqrt(3/4) * self.g3 * self.a * k

        # g4
        H13 = H31 = H42 = H24 = - (3/4)**(1/2) * self.a * self.g4 * k

        H = np.array([  [H11, H12, H13, H14],
                        [H21, H22, H23, H24],
                        [H31, H32, H33, H34],
                        [H41, H42, H43,H44]]).squeeze()
        return H

    ######################
    ### Band Structure ###
    ######################

    def Dispersion(self,k,u,band,approx='Common'):
        '''
        Returns the energy (J) of an electron with wavevector k (rad/m)
        in first (band=1) or second (band=2) conduction band.
        Only approximation is g3=0.
        To get valence bands, simply result multiply by -1.
        '''
        p = sc.hbar * self.vF * k

        if approx == 'Common':
            radical=(self.g1**4)/4 + (u**2 + self.g1**2)*(p**2)
            return np.sqrt( (self.g1**2)/2 + (u**2)/4 + p**2 + ((-1)**(band))*np.sqrt(radical) )

        if approx == 'LowEnergy':
            '''
            Low Energy effective. Eigenvalues of 

            H = ( ( u/2, p^2 / 2m ) , ( p^2/2m, -u/2 ) )

            '''

            meff = ( self.g1 / (2 * (self.vF)**2) )
            return np.sqrt( (u/2)**2 + ( (hbar * k)**2 / (2*meff) )**2 )

        if approx == 'None':
            '''
            No approximation. Compute eigenvalues of Hamiltonian
            '''
            k = k.squeeze()
            u = u.squeeze()
            disp = np.empty(np.shape(k))

            for i, wn in enumerate(k):
                disp[i] = linalg.eigvalsh(self.Hamiltonian(wn,u))[1+band]

            return np.array(disp).squeeze()

    def kmin(self,u, band=1):
        '''
        Returns positive wavenumber at the minimum of the first band in 1/m.

        Parameters
        ----------
        u :     array-like
                Interlayer potential energy difference in units J.

        band:   First (second) conduction band 1 (2).
        '''
        k2 = ( u**2 / (2*hbar*self.vF)**2 ) * ( (2*self.g1**2 + u**2) /( self.g1**2 + u**2 ) )
        return np.sqrt(k2)

    def emin(self,u):
        '''
        Returns minimum of the first band in Joules.
        '''
        emin2 = (u/2)**2 * (self.g1**2 / ( self.g1**2 + u**2 ) )
        return np.sqrt(emin2)

    def DOS(self, e, u):
        '''
        Returns the density of states per unit area (1/m^2) as 
        a function of energy given the gap u
        '''
        e = np.atleast_1d(abs(e))
        
        # Define the multiplicative factor out front
        # Set to 0 is energy is below the minimum
        mult = (e>self.emin(u)) * ( e / (pi * hbar**2 * self.vF**2) )
        
        # Calculate the discriminant
        # Remember, we wil need to divide by it's sqrt
        # So set values disc<=0 to 1
        # We will multiply the result by zero for these energies anyway later on.
        disc = e**2 * (self.g1**2 + u**2) - self.g1**2 * u**2 / 4
        disc = (e>self.emin(u))*disc + (e<=self.emin(u))*1
        
        # Calculate quantities proportional to derivatives of k^2
        propdkp2 = 2 + (self.g1**2 + u**2)/np.sqrt(disc)
        propdkm2 = 2 - (self.g1**2 + u**2)/np.sqrt(disc)
        
        # If energy is above sombrero region, add the positive solution
        # If within, also subtract the negative solution
        propdos = (e>self.emin(u))*propdkp2 - (e<=abs(u/2))*propdkm2
        return (mult * propdos)

    def Pdiff(self,k,vminus,approx='Common'):
        '''Returns the probability difference between finding an ELECTRON on the TOP layer minus the BOTTOM layer.'''
        
        u = 2*q*(vminus+np.sign(vminus)*0.0000001)
        
        if approx=='Common':
            e = self.Dispersion(k,u,1)

            K = hbar*self.vF*(k+1)
            
            numerator = (e**2 - u**2/4)**2 + 4*K**2*e**2 - K**4
            denominator = (e**2 - u**2/4)**2 + K**2*u**2 - K**4
            
            return - ( u / (2*e) ) * ( numerator / denominator )

        if approx=='LowEnergy':
            meff = ( self.g1 / (2 * (self.vF)**2) )
            denominator_squared = ( ( (hbar*k)**2/meff )**2 + u**2 )
            
            
            return - u / np.sqrt(denominator_squared)

        if approx=='None':
            k = np.atleast_1d(k).squeeze()
            u = np.atleast_1d(u).squeeze()
            deltapsi = []
            # Eigenvectors of 
            for i,wn in enumerate(k):
                v = linalg.eigh( self.Hamiltonian(wn,u) )[1]

                psi = v[:,-2] # Second highest band (first conduction)

                deltapsi.append(psi[0]**2 + psi[1]**2 - psi[2]**2 - psi[3]**2)

            return np.array(deltapsi).squeeze()

    def kFermi(self,n,u,pm):
        '''
        Returns Fermi vector kF+ for pm=1 and kF- for pm=2 in units rad/m
        '''
            
        # Define the more complicated factors and terms
        numerator = (pi * hbar**2 *self.vF**2 * n)**2 + ( self.g1*u )**2
        denominator = self.g1**2 + u**2
        pmterm = 2*pi*hbar**2 * self.vF**2 * abs(n) # abs for fact that electrons and holes symmetric
        
        # Factor proportional to k**2
        propk2 = ( numerator / denominator ) + u**2 + (-1)**(pm-1) * pmterm
        
        # If the fermi level is above u/2, set kF- to zero
        # This says that the region of occupied states is now a disk
        if pm%2==0:
            propk2 = (propk2 >= 0) * propk2
            propk2 = (self.Dispersion(self.kFermi(n,u,1),u,1)<u/2) * propk2
        
        return np.sqrt( propk2 ) / (2*hbar*self.vF)

    def eFermi(self,n,u):
        '''
        Returns the Fermi level (Joules) given density n and interlayer potential energy difference u
        Positive n returns a positive Fermi level, meaning positive carrier densities are electrons by convention.
        '''
        
        numerator = (hbar**2 * self.vF**2 * n *pi)**2 + (self.g1 * u)**2
        denominator = 4 * (self.g1**2 + u**2)
        
        return np.sign(n) * np.sqrt( numerator / denominator )

    #########################
    ### Carrier Densities ###
    #########################

    def nplusT0(self,vplus,vminus,approx='Fermi'):
        """
        Analytically computes the electron density at zero temperature.
        Faster than Bilayer.nplus() since this function allows
        for vectorized operations.
        """

        # Convert voltages to energies
        eF = eVtoJ*vplus
        u  = 2*eVtoJ*vminus

        if approx == 'Fermi':
            term1 = 4*(self.g1**2 + u**2)*(eF**2)
            term2 = -(self.g1**2)*(u**2)

            prop = (hbar**2 * self.vF**2 * pi)**(-1)

            n = np.sign(eF)*prop*np.sqrt( term1 + term2 )

        if approx == 'Common':
            # Calculate the radical
            radical = (self.g1**2+u**2) * eF**2 - self.g1**2 * u**2 / 4

            # For energies within the gap, radical is negative, so set it to 0 instead
            radical = (radical>=0)*radical

            # Proportional to the square of the Fermi wavevectors
            kFp2 = (eF**2 + u**2/4) + np.sqrt(radical)
            kFm2 = (eF**2 + u**2/4) - np.sqrt(radical)

            # For kFm2, if eF > u/2, set to zero
            kFm2 = (abs(eF) <= abs(u/2)) * kFm2
            
            # Calculate the proportionality factor
            # Includes:
            #     1/(hbar vF)**2 from formula for kF
            #     1/pi from n = (kFp2 - kFm2)/pi
            #     Sets to zero if Fermi in the gap
            prop = (abs(eF)>self.emin(u))*np.sign(eF)*(1 / (hbar**2 * self.vF**2 * pi))

            n = prop * (kFp2 - kFm2)

            return n

        if approx == 'LowEnergy':
            """
            See Young and Levitov 2011.
            """
            meff = ( self.g1 / (2 * (self.vF)**2) )

            nu0 = 2 * meff * q  / (pi * hbar**2)

            energy_diff = (np.abs(eF)>np.abs(u/2)) * (eF**2 - (u/2)**2)
            return (nu0/q) * np.sign(eF) * np.sqrt(energy_diff)

    def nminusT0(self,vplus,vminus):

        meff = ( self.g1 / (2 * (self.vF)**2) )
        nu0 = 2 * meff * q  / (pi * hbar**2)

        prop = nu0 * vminus
        
        # Find cutoff energy. Approximate it as the vminus=0 case
        Lambda = self.Dispersion( 1 / (np.sqrt(3) * self.a), -2*q*vminus, 1 ) / q
        
        # Compute the denominator of the log
        metal = abs(vplus) >= abs(vminus)
        den = (metal) * np.abs(vplus) + np.sqrt( metal * vplus**2 + (-1)**metal * vminus**2 ) 
        
        return prop * np.log(2 * Lambda / den)


    #################
    ### Screening ###
    #################



    def screened_vminus2(self,nplus,vminus):
        """
        The screened value of vminus given the total charge nplus
        """
        a, b = -1, 1

        vminus_screened = []

        for vm in vminus:
            vm0 = vm
            vp0 = self.eFermi(nplus, -2*q*vm) / q

            def f1(vm1):
                return (vm1 - vm) + (q / (4*self.C))*self.nminus(vp0,vm1,0)

            vm1 = optimize.brentq(f1,a,b)
            vp1 = self.eFermi(nplus, -2*q*vm1) / q

            while (vm1-vm0)**2 + (vp1-vp0)**2 > 0.0001:
                vp0, vm0 = vp1, vm1

                def f1(vm1):
                    return (vm1 - vm) + (q / (4*self.C))*self.nminus(vp0,vm1,0)

                vm1 = optimize.brentq(f1,a,b)
                vp1 = self.eFermi(nplus, -2*q*vm1) / q
            
            vminus_screened.append(vm1)

        return np.array(vminus_screened)

    def screened_newton(self,vplus,vminus):
        n = self.nplus(vplus,vminus,0)

        def f1(v):
            return (v[1] - vminus) + (q / (4*self.C))*self.nminusT0(v[0],v[1])

        def f2(v):
            return n - self.nplus(v[0],v[1],0)

        v = Newton.Newton2D(f1,f2,np.array([vplus,vminus]))

        return v


    ##########################
    ### Plotting Functions ###
    ##########################

    def plot_band_structure(self,n=None,vplus=None,vminus=None,schematic=False,savefile=None):
        '''

        Parameters
        ----------

        n:      Scalar; the charge area density in m^-2

        vplus:  Scalar; the potential of the BLG. Related to fermi level by vplus = eF/|e|

        vminus: Scalar; half the potential difference between the layers.

        schematic: Bool; Setting to True makes it look more like a schematic than a graph.

        savefile;  str; 
        '''


        if n and vplus:
            print("Can't specify n and vplus simultaneously")
            return

        u = 2*q*vminus

        if n or n==0:
            vplus = self.eFermi(n,u)/q
        if vplus or vplus==0:
            n = self.nplusT0(vplus,vminus)


        kmax = 5e8
        k = np.linspace(-kmax,kmax,num=100)

        en_con = (self.Dispersion(k,u,band=1)/q)
        en_val = (-self.Dispersion(k,u,band=1)/q)
        fig, ax = plt.subplots()

        # Plots to ensure proper plotting
        ax.plot(k*1e-10,200*np.ones_like(k),'w-')
        ax.plot(k*1e-10,-200*np.ones_like(k),'w-')
        ax.set_ylim(-201,201)
        ax.plot(k*1e-10,en_con*1e3,'k-')
        ax.plot(k*1e-10,en_val*1e3,'k-')
        #ax.plot(k*1e-10,np.zeros_like(k),'k--',label='$V_+$')


        if n>0:
            #Fill the entire lower band
            ax.fill_between(k*1e-10,en_val[0]*1e3*np.ones_like(en_val),en_val*1e3,where=en_val>en_val[0],facecolor='b')
            # Then fill portion of upper band
            ax.fill_between(k*1e-10,en_con*1e3,vplus*1e3*np.ones_like(k),where=vplus>en_con,facecolor='b')
        if n<0:
            # Fill entire lower band
            ax.fill_between(k*1e-10,en_val[0]*1e3*np.ones_like(en_val),en_val*1e3,where=en_val>en_val[0],facecolor='b')
            # then fill upper portion of lower band with white
            ax.fill_between(k*1e-10,en_val*1e3,vplus*1e3*np.ones_like(k),where=vplus<en_val,facecolor='w')            

        ax.set_xlabel('k (1/A)')
        ax.set_ylabel('Energy (meV)')

        if schematic:
            ax.set_axis_off()

        if savefile:
            plt.savefig(savefile,dpi=150,bbox_inches='tight');

        plt.show()


    ##################
    ### OLD Screening METHOD ###
    ##################

    def nplus(self,vplus,vminus, T, approx='Common',points = 10000):
        '''
        Returns the electron carrier density for various electrostatic potentials vplus, vminus.
        Convention is that electrons have positive carrier density while holes have negative.
        '''

        # Treat inputs as ndarrays so we can take advantage of broadcasting
        vplus = np.atleast_1d(vplus)
        vminus = np.atleast_1d(vminus)

        vplus = vplus.reshape(1,1,len(vplus))
        vminus = vminus.reshape(1,len(vminus),1)

        # Domain over first Brillouin zone
        ks = np.linspace(0,1/(np.sqrt(3)*self.a), num=points).reshape((points,1,1))

        # Calculate the kinetic energy
        KE = self.Dispersion(ks, -2*q*vminus,1,approx)

        # Evaluate Fermi-Dirac
        FD = (sd.FermiDirac(KE-q*vplus,T)-sd.FermiDirac(KE+q*vplus,T))

        # Define integrand
        integrand = ( 2 / np.pi ) * ks * FD

        return np.squeeze(integrate.trapz(integrand,ks,axis=0))

    def nminus(self,vplus,vminus, T, approx='Common', points=10000):
        '''
        Returns the electron carrier density for various electrostatic potentials vplus.
        Convention is that electrons have positive carrier density while holes have negative.
        '''

        if approx == 'None':
            print('Not yet supported')
            return
        # Treat inputs as ndarrays so we can take advantage of broadcasting
        vplus = np.atleast_1d(vplus)
        vminus = np.atleast_1d(vminus)

        vplus = vplus.reshape(1,1,len(vplus))
        vminus = vminus.reshape(1,len(vminus),1)

        # Domain over first Brillouin zone
        ks = np.linspace(0,1/(np.sqrt(3)*self.a), num=points).reshape((points,1,1))

        # Calculate the kinetic energy
        KE = self.Dispersion(ks, -2*q*vminus,1, approx)

        # Evaluate Fermi-Dirac
        # Minus sign comes from...
        FD = (sd.FermiDirac(KE-q*abs(vplus),T))#-Temperature.FermiDirac(-KE-q*vplus,T)

        # Define integrand
        integrand =  ( 2 /np.pi ) * ks * self.Pdiff(ks,vminus,approx='LowEnergy') * FD

        nm = np.squeeze(integrate.trapz(integrand,ks,axis=0))

        return nm

    def generate_nplus_nminus(self,vplus,vminus,T):
        """
        Generates and saves high-resolution surfaces of nplus(vplus,vminus)
        and nminus(vplus,vminus). Only generate for first quadrant (vplus,vminus > 0)
        since surfaces have symmetry properties.
        """
        save_dir = os.path.join(self.this_dir,
                                'CarrierDensities',
                                'Temp_{:.2E}'.format(T))

        if os.path.exists(save_dir):
            print('Carrier densities for T = {} K have already been generated'.format(T))
            
        #else:
        #    os.makedirs(save_dir)

        if np.any(vplus< 0) or np.any(vminus<0):
            print('Some voltages were negative in the ranges\n')
            print(  vplus[0].squeeze(),
                    ' < vplus < ',
                    vplus[-1].squeeze(),
                    ' {} points'.format(np.shape(vplus)[0]))
            print(  vminus[0].squeeze(),
                    ' < vminus < ', vminus[-1].squeeze(),
                    ' {} points'.format(np.shape(vminus)[0]))

            vplus   = np.linspace(0,vplus[-1],num=np.shape(vplus)[0])
            vminus  = np.linspace(0,vminus[-1],num=np.shape(vminus)[0]).reshape(np.shape(vminus))

            print('\nInstead, generating over the ranges\n')
            print(  '0 < vplus < ', vplus[-1].squeeze(),
                    ' {} points'.format(np.shape(vplus)[0]))
            print(  '0 < vminus< ', vminus[-1].squeeze(),
                    ' {} points'.format(np.shape(vminus)[0]))
            print()

        # Choose the size of the batches we will generate
        d = 10

        # Check that it is compatible with the lengths of vplus and vminus
        if len(vplus) % d != 0 or len(vminus) % d != 0:
            print('Batch size (d) incompatible with voltage arrays')
            print('d= {} does not evenly divide either len(vplus)= {} or len(vminus) = {}'.format(d,len(vplus),len(vminus)))
            return

        nplus_surface = np.empty(np.shape(vminus*vplus))
        nminus_surface = np.empty(np.shape(vminus*vplus))

        for i in range(int(len(vplus)/d)):
            i_frac = i / int(len(vplus)/d)
            for j in range(int(len(vminus)/d)):
                j_frac  = (j / int(len(vminus)/d)) * (1 / int(len(vplus)/d))
                percentage = round(100* (i_frac + j_frac),2)
                print('{} % Finished'.format(percentage))
                nplus_surface[d*j:d*j+d,d*i:d*i+d]=self.nplus(vplus[d*i:d*i+d],vminus[d*j:d*j+d,:],T)
                nminus_surface[d*j:d*j+d,d*i:d*i+d]=self.nminus(vplus[d*i:d*i+d],vminus[d*j:d*j+d,:],T)

        # Save the surfaces
        np.save(save_dir+'nplus_surface.npy',nplus_surface)
        np.save(save_dir+'nminus_surface.npy',nminus_surface)

        # Save the voltages
        np.save(save_dir+'vplus.npy', vplus)
        np.save(save_dir+'vminus.npy', vminus)

    def get_vplus(self,T):
        """
        Returns the vplus array saved in ...
        Doubles the range to negative values
        """
        save_dir = os.path.join(self.this_dir,
                                'CarrierDensities',
                                'Temp_{:.2E}'.format(T))
        vplus = np.load(save_dir+'vplus.npy')
        return np.concatenate((-vplus[:0:-1],vplus))

    def get_vminus(self,T):
        save_dir = os.path.join(self.this_dir,
                                'CarrierDensities',
                                'Temp_{:.2E}'.format(T))
        vminus = np.load(save_dir+'vminus.npy')
        return np.concatenate((-vminus[:0:-1],vminus))

    def get_nplus(self,T):
        save_dir = os.path.join(self.this_dir,
                                'CarrierDensities',
                                'Temp_{:.2E}'.format(T))
        nplus_surface = np.load(save_dir+'nplus_surface.npy')
        nplus_surface = np.concatenate((nplus_surface[:0:-1,:],nplus_surface))
        nplus_surface = np.concatenate((-nplus_surface[:,:0:-1],nplus_surface),axis = 1)
        return nplus_surface

    def get_nminus(self,T):
        save_dir = os.path.join(self.this_dir,
                                'CarrierDensities',
                                'Temp_{:.2E}'.format(T))
        nminus_surface = np.load(save_dir+'nminus_surface.npy')
        nminus_surface = np.concatenate((-nminus_surface[:0:-1,:],nminus_surface))
        nminus_surface = np.concatenate((nminus_surface[:,:0:-1],nminus_surface),axis = 1)
        return nminus_surface

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
