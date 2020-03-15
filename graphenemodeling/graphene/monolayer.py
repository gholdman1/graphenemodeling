"""
======================================================
Monolayer (:mod:`graphenemodeling.graphene.monolayer`)
======================================================

Functions
=========

Band structure
--------------

.. toctree::
    :maxdepth: 1

    graphene.monolayer.Hamiltonian
    graphene.monolayer.CarrierDispersion
    graphene.monolayer.DensityOfStates
    graphene.monolayer.FermiWavenumber
    graphene.monolayer.CarrierDensity
    graphene.monolayer.ChemicalPotential

Optical Properties
------------------

.. toctree::
    :maxdepth: 1

    graphene.monolayer.Polarizibility
    graphene.monolayer.OpticalConductivity


Plasmonics
----------

.. toctree::
    :maxdepth: 1

    graphene.monolayer.PlasmonDispersion

Models
======

Different models for graphene are valid in different regimes.
Many of the functions provided here have a ``model`` option allowing the use to choose which to use.

Full Tight Binding (:mod:`model='FullTightBinding'`)
----------------------------------------------------
The tight binding model Hamiltonian is

    .. math::

        H = \\gamma_0 \\left(\\array{
                                0 & f(k) \n
                                f(k)^* & 0
                                        } \\right)

where :math:`f(k)= e^{ik_x a/2} + 2 e^{-i k_x a/ 2}\\cos(k_y a \\sqrt{3}/2)`


Low Energy (:mod:`model='LowEnergy'`)
-------------------------------------
An low-energy approximation to the Tight Binding model is what we are calling the Low Energy model.
This is valid for small momenta, i.e. :math:`\\hbar v_F k\\ll \\gamma_0`.


References
==========

[1] Castro Neto, A.H., Guinea, F., Peres, N.M.R., Novoselov, K.S., and Geim, A.K. (2009).
The electronic properties of graphene. Rev. Mod. Phys. 81, 109–162.
https://link.aps.org/doi/10.1103/RevModPhys.81.109

"""

import numpy as np
import scipy.constants as sc
from scipy import special, optimize, integrate

from graphenemodeling.graphene.base import BaseGraphene
import graphenemodeling.graphene._constants as _c

import graphenemodeling.statistical_distributions as sd

############
# Geometry #
############

def UnitCell(m,n):
    '''Positions of unit cell

    Parameters
    ----------

    m, n:   Unit cell indices.

    References
    ----------

    [1] Castro Neto, A.H., Guinea, F., Peres, N.M.R., Novoselov, K.S., and Geim, A.K. (2009). The electronic properties of graphene. Rev. Mod. Phys. 81, 109–162. https://link.aps.org/doi/10.1103/RevModPhys.81.109.

    '''
    a1 = np.array(_c.a1)
    a2 = np.array(_c.a2)

    pos = m*a1 + n*a2

    return pos

def AtomicPosition(m,n,i):

    delta = _c.a/2*np.array([1,3**(1/2)])
    pos = UnitCell(m,n)+i*delta

    return pos

##################
# Band Structure #
##################

def Hamiltonian(k,model='LowEnergy'):
    '''Hamiltonian in momentum space.

    Parameters
    ----------

    k:      array-like, wavevector of carrier. Use complex k=kx + iky for 2D wavevectors.

    Returns
    ----------

    H:      2x2 array, Tight-binding Hamiltonian evaluated at k.

    Notes
    -----

    Let :math:`k=k_x+ik_y`. Then the common ``model=LowEnergy`` approximation is

    .. math::

        H = \\hbar v_F\\left(\\array{
                                    0 & k \n
                                    k^* & 0

                                    } \\right)

    while the ``model=FullTightBinding`` expression is given by 

    .. math::

        H = \\gamma_0 \\left(\\array{
                                0 & f(k) \n
                                f(k)^* & 0
                                        } \\right)

    where :math:`f(k)= e^{ik_x a/2} + 2 e^{-i k_x a/ 2}\\cos(k_y a \\sqrt{3}/2)`

    References
    ----------

    [1] Slonczewski, J.C., and Weiss, P.R. (1958). Band Structure of Graphite. Phys. Rev. 109, 272–279. https://link.aps.org/doi/10.1103/PhysRev.109.272.

    [2] Falkovsky, L.A., and Varlamov, A.A. (2007). Space-time dispersion of graphene conductivity. Eur. Phys. J. B 56, 281–284. https://link.springer.com/article/10.1140/epjb/e2007-00142-3.

    [3] Christensen, T. (2017). From Classical to Quantum Plasmonics in Three and Two Dimensions (Cham: Springer International Publishing). http://link.springer.com/10.1007/978-3-319-48562-1.



    '''

    if model == 'LowEnergy':
        H11 = 0
        H12 = sc.hbar * _c.vF * k
        H21 = np.conj(H12)
        H22 = 0

    if model == 'FullTightBinding':
        kx = np.real(k)
        ky = np.imag(k)

        H11 = 0
        H12 = _c.g0 * (   np.exp(1j*kx*_c.a/2)
                            + 2*np.exp(-1j*kx*_c.a/2)*np.cos(ky*_c.a*np.sqrt(3)/2) )
        H21 = np.conj(H12)
        H22 = 0

    H = np.array( [[H11, H12],
                    [H12, H22] ])

    return H

def CarrierDispersion(k,model,eh=1,g0prime=_c.g0prime):
    '''The dispersion of Dirac fermions in monolayer graphene.

    These are the eigenvalues of the Hamiltonian.
    However, in both the ``LowEnergy`` model and the ``FullTightBinding`` model, we use closed form solutions rather than solving for the eigenvalues directly.
    This saves time and make broadcasting easier.

    Parameters
    ----------

    k:          array-like, complex
                Wavevector of Dirac fermion relative to K vector.
                For 2D wavevectors, use :math:`k= k_x + i k_y`.

    model:      string
                ``'LowEnergy'``: Linear approximation of dispersion.

                ``'FullTightBinding'``: Eigenvalues of tight-binding approximation. We use a closed form rather than finding eigenvalues of Hamiltonian to save time and avoid broadcasting issues.

    eh:         int
                Band index: 
                ``eh=1``  returns conduction band,
                ``eh=-1`` returns valence band

    Returns
    ----------

    array-like

    Notes
    -----

    When ``model='LowEnergy'``,

    .. math::

        E =\\pm\\hbar v_F |k|

    When ``model=FullTightBinding``,

    .. math::

        E = \\pm \\gamma_0 \\sqrt{3 + f(k)} - \\gamma_0'f(k)

    where :math:`f(k)= 2 \\cos(\\sqrt{3}k_y a) + 4 \\cos(\\sqrt{3}k_y a/2)\\cos(3k_xa/2)`.

    Both expressions are equivalent to diagonalizing the Hamiltonian of the corresponding ``model``.

    Examples
    --------

    Plot the Fermion dispersion relation.

    .. plot::

        >>> import matplotlib.pyplot as plt
        >>> from graphenemodeling.graphene import monolayer as mlg
        >>> from scipy.constants import elementary_charge as eV
        >>> eF = 0.4*eV
        >>> kF = mlg.FermiWavenumber(eF,model='LowEnergy')
        >>> k = np.linspace(-2*kF,2*kF,num=100)
        >>> conduction_band = mlg.CarrierDispersion(k,model='LowEnergy')
        >>> valence_band = mlg.CarrierDispersion(k,model='LowEnergy',eh=-1)
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

    Plot the full multi-dimensional dispersion relation with a particle-hole asymmetry. Replicates Figure 3 in Ref. [1].

    .. plot::

        >>> from graphenemodeling.graphene import monolayer as mlg
        >>> from graphenemodeling.graphene import _constants as _c
        >>> import matplotlib.pyplot as plt
        >>> from mpl_toolkits import mplot3d # 3D plotting
        >>> kmax = np.abs(_c.K)
        >>> emax = mlg.CarrierDispersion(0,model='FullTightBinding',g0prime=-0.2*_c.g0)
        >>> kx = np.linspace(-kmax,kmax,num=100)
        >>> ky = np.copy(kx)
        >>> k = (kx + 1j*ky[:,np.newaxis]) + _c.K # k is relative to K. Add K to move to center of Brillouin zone
        >>> conduction_band = mlg.CarrierDispersion(k,model='FullTightBinding',eh=1,g0prime=-0.2*_c.g0)
        >>> valence_band = mlg.CarrierDispersion(k,model='FullTightBinding',eh=-1,g0prime=-0.2*_c.g0)
        >>> fig = plt.figure(figsize=(8,8))
        >>> fullax = plt.axes(projection='3d')
        >>> fullax.view_init(20,35)
        >>> KX, KY = np.meshgrid(kx,ky)
        >>> fullax.plot_surface(KX/kmax,KY/kmax,conduction_band/_c.g0,rstride=1,cstride=1,cmap='viridis',edgecolor='none')
        <...
        >>> fullax.plot_surface(KX/kmax,KY/kmax,valence_band/_c.g0,rstride=1,cstride=1,cmap='viridis',edgecolor='none')
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

    References
    ----------

    [1] Castro Neto, A.H., Guinea, F., Peres, N.M.R., Novoselov, K.S., and Geim, A.K. (2009).
    The electronic properties of graphene. Rev. Mod. Phys. 81, 109–162. 
    https://link.aps.org/doi/10.1103/RevModPhys.81.109.


    '''

    if eh!=1 and eh!=-1:
        raise ValueError('eh must be either 1 or -1')

    if model == 'LowEnergy':

        return eh*sc.hbar*_c.vF*np.abs(k)

    if model == 'FullTightBinding':

        k = k - _c.K
        f = lambda k: (2*np.cos(np.sqrt(3)*np.imag(k)*_c.a) 
                        + 4*np.cos((np.sqrt(3)*np.imag(k)/2)*_c.a)*np.cos((3/2)*np.real(k)*_c.a) )

        # [sic] eh only applies to first term
        return eh*_c.g0*np.sqrt(3+ f(k)) - g0prime*f(k)

def FermiWavenumber(eFermi,model,g0prime=_c.g0prime):
    '''
    The Fermi wavenumber, i.e. the wavenumber of the state at
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
    >>> kF = mlg.FermiWavenumber(eF, model='LowEnergy')
    >>> mlg.CarrierDispersion(kF,model='LowEnergy')/eV
    0.4
    '''

    if model == 'LowEnergy':
        return np.abs(eFermi) / (sc.hbar*_c.vF)

    if model == 'FullTightBinding':
        '''
        Need to finish off this code-finding procedure
        '''
        eh = np.sign(eFermi)

        # f is zero when kf is correct value
        f = lambda kf: eFermi - CarrierDispersion(kf, model='FullTightBinding',eh=eh,g0prime=g0prime)

        # Choose LowEnergy answer for initial starting point
        kf0 = FermiWavenumber(eFermi,model='LowEnergy',g0prime=g0prime)

        result = optimize.root_scalar(f,x0=kf0,x1=kf0*.9,rtol=1e-10).root

        return result

def DensityOfStates(E,model,g0prime=_c.g0prime):
    '''
    The density of states per square meter of graphene at energy :math:`E`.

    Parameters
    ----------

    E:  array-like
        Energy :math:`E` at which to evaluate density of states.

    model:      string
                ``'LowEnergy'``: Linear approximation of dispersion.

                ``'FullTightBinding'``: Eigenvalues of tight-binding approximation. We use a closed form rather than finding eigenvalues of Hamiltonian to save time and avoid broadcasting issues.

    Returns
    -------

    array-like
        Density of states, units are states per J-m^2

    Notes
    -----

    For ``model==LowEnergy``, the form is simply

    .. math::

        \\rho(E)=\\frac{2}{\\pi}\\frac{|E|}{\\hbar^2 v_F^2}

    whereas the ``FullTightBinding`` model has a much more complicated form (eqn. 14 of [2])

    .. math::

        \\rho(E)=\\frac{4}{\\pi^2}\\frac{|E|}{\\gamma_0^2}\\frac{1}{\\sqrt{Z_0}}\\mathbf{F}\\left(\\frac{\\pi}{2},\\sqrt{\\frac{Z_1}{Z_0}}\\right)

    where :math:`\\mathbf{F}(\\pi/2,x)` is the complete elliptic integral of the first kind (see `scipy.special.ellipk`_) and

    .. _scipy.special.ellipk: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipk.html

    .. math::

        Z_0 = \\left\{\\array{    
                    (1 + |E/\\gamma_0|)^2 - \\frac{[(E/\\gamma_0)^2-1]^2}{4}, & |E|\\leq \\gamma_0 \n
                    4|E/\\gamma_0|, & -3\\gamma_0\\leq E \\leq -\\gamma_0,\\gamma_0\\leq E\\leq 3\\gamma_0
                }\\right.

        Z_1 = \\left\{\\array{
                    4|E/\\gamma_0|, & |E|\\leq \\gamma_0 \n    
                    (1 + |E/\\gamma_0|)^2 - \\frac{[(E/\\gamma_0)^2-1]^2}{4}, & -3\\gamma_0\\leq E \\leq -\\gamma_0,\\gamma_0\\leq E\\leq 3\\gamma_0
                }\\right.

    Examples
    --------
    Plot the density of states for ``model=LowEnergy`` approximation and ``model=FullTightBinding`` model. Replicates Fig. 5 in Ref. [2].

    .. plot::

        >>> from graphenemodeling.graphene import monolayer as mlg
        >>> from graphenemodeling.graphene import _constants as _c
        >>> import matplotlib.pyplot as plt
        >>> E = np.linspace(-3,3,num=200) * _c.g0
        >>> DOS_low = mlg.DensityOfStates(E,model='LowEnergy')
        >>> DOS_full = mlg.DensityOfStates(E,model='FullTightBinding')
        >>> plt.plot(E/_c.g0,DOS_full/np.max(DOS_full),'k-',label='FullTightBinding')
        [<...
        >>> plt.plot(E/_c.g0,DOS_low/np.max(DOS_full),'k-.',label='LowEnergy')
        [<...
        >>> plt.xlabel('$E/\\gamma_0$')
        Text...
        >>> plt.ylabel('DOS (a.u.)')
        Text...
        >>> plt.legend()
        <...
        >>> plt.show()
    
    References
    ----------

    [1] Hobson, J.P., and Nierenberg, W.A. (1953). The Statistics of a Two-Dimensional, Hexagonal Net. Phys. Rev. 89, 662–662. https://link.aps.org/doi/10.1103/PhysRev.89.662.

    [2] Castro Neto, A.H., Guinea, F., Peres, N.M.R., Novoselov, K.S., and Geim, A.K. (2009).
    The electronic properties of graphene. Rev. Mod. Phys. 81, 109–162.
    https://link.aps.org/doi/10.1103/RevModPhys.81.109.

    '''

    if model=='LowEnergy':

        E = np.abs(E)

        DOS = 2 * E / (sc.pi*(sc.hbar*_c.vF)**2)

        return DOS

    elif model=='FullTightBinding':

        if g0prime!=0:
            raise Exception('Not supported for g0prime!=0.\nSetting g0prime=0')
            g0prime=0


        prefactor = 4*np.abs(E) / (sc.pi*_c.g0)**2

        def fZ0(E):
            if np.abs(E)<np.abs(_c.g0):
                term1 = (1+np.abs(E/_c.g0))**2
                term2 = -((E/_c.g0)**2 - 1)**2 / 4

                return term1 + term2

            else: return 4*np.abs(E/_c.g0)

        def fZ1(E):
            if np.abs(E)<np.abs(_c.g0):
                return 4*np.abs(E/_c.g0)

            else:
                term1 = (1+np.abs(E/_c.g0))**2
                term2 = -((E/_c.g0)**2 - 1)**2 / 4

                return term1 + term2

        dos = np.empty_like(E)

        for j, e in np.ndenumerate(E):
            Z0 = fZ0(e)
            Z1 = fZ1(e)
            ellip = special.ellipk(np.sqrt(Z1/Z0))

            dos[j] = (1/np.sqrt(Z0)) * ellip
        DOS = prefactor * dos /_c.A

        return DOS
    else:
        print('The model %s is not available' % (model))

def CarrierDensity(mu,T,model):
    '''
    Computes the carrier density directly from the band structure.

    Parameters
    ----------

    mu: array-like
        Chemical potential :math:`\\mu`

    T:  scalar
        Temperature :math:`T`

    Returns
    -------

    array-like

    Notes
    -----

    When ``T>0``, we use the standard formula of integrating the Fermi-Dirac distribution over the density of states :math:`\\rho`  .

    .. math::

        n = \\int_{-\\infty}^{\\infty} \\frac{1}{e^{(\\epsilon-\\mu)/k_BT}+1}\\rho(\\epsilon) d\\epsilon

    For graphene at ``T=0``, this reduces to

    .. math::

        n = \\frac{\\mu^2}{\\pi\\hbar^2 v_F^2}

    '''

    if T<0:
        raise ValueError('Temperature T must be nonnegative')
    if T==0 and model=='LowEnergy':
        eFermi=mu # chemical potential at T=0 is called Fermi level
        n = (eFermi / (sc.hbar*_c.vF))**2 / np.pi
        return n

    if T>0:
        n = np.empty_like(mu)
        for i, m in np.ndenumerate(mu):
            p_electron = lambda e: DensityOfStates(e,model) * sd.FermiDirac(e-m,T)
            p_hole = lambda e: DensityOfStates(e,model) * (1 - sd.FermiDirac(e-m,T))
            n[i] = ( integrate.quad(p_electron,0,3*_c.g0,points=(_c.g0,m))[0] -
                     integrate.quad(p_hole,-3*_c.g0,0,points=(-_c.g0,-m))[0]   )
    return n

def ChemicalPotential(n,T=0,model='LowEnergy'):
    '''Returns the chemical potential given the carrier density.

    Essentially the inverse of Carrier Density.

    Parameters
    ----------

    n:  array-like
        Carrier density in units m^-2.

    T:  scalar
        Temperature (K)

    Returns
    -------

    array-like

    Notes
    -----

    When ``T=0`` and ``model='LowEnergy'`` simultaneously, a closed form expression is used.

    .. math::

        E_F = \\hbar v_F\\sqrt{\\pi n}
    
    For ``T>0``, we use a numerical routine regardless of model.

    '''

    if T==0 and model=='LowEnergy':
        return sc.hbar*_c.vF*np.sqrt(sc.pi*n)

    else:
        ## Numerically solve

        # f is zero when n is correct
        f = lambda mu: n - CarrierDensity(mu,T,model=model)

        # Use T=0 value to estimate start
        mu0 = ChemicalPotential(n,T=0,model=model)

        result = optimize.root_scalar(f,x0=mu0,x1=mu0*1.1+0.1*sc.elementary_charge,
                                        rtol=1e-10).root

        return result

########################
# Electrical Transport #
########################

def Mobility(Te,mu0,T0):
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

def ScatteringRate(mobility,eFermi):
    '''
    Estimated DC scattering rate from mobility.

    Parameters
    ----------

    mobility:   scalar, mobility (m^2/V-s)

    eFermi:     scalar, Fermi level (J)

    Returns
    ----------

    rate:       scalar, scattering rate
    '''

    # Scattering time
    tau = mobility*eFermi / (sc.elementary_charge*_c.vF**2)

    rate = 1/tau
    return rate

def Polarizibility(q,omega,gamma,eFermi,T=0):
    '''
    The Polarizibility :math:`\\chi^0`` of graphene.

    Parameters
    ----------

    q:      array-like, difference between scattered wavevector and incident

    omega:  array-like, frequency

    gamma:  scalar, scattering rate due to mechanisms such as impurities (i.e. not Landau Damping)
                    We use the Mermin-corrected Relaxation time approximation (Eqn 4.9 of Ref 1)

    eFermi: scalar, Fermi level of graphene.

    Notes
    -----
    The Polarizibiliy function in graphene. Can be derived from a
    self-consistent field method or the Kubo formula.

    .. math::

        \\chi^0(\\mathbf q, \\omega) = \\frac{g}{A}\\sum_{nn'\\mathbf k}\\frac{f_{n\\mathbf k} - f_{n'\\mathbf{k+q}}}{\\epsilon_{n\\mathbf k}-\\epsilon_{n'\\mathbf{k+q}}+\\hbar(\\omega+i\\eta)}|\\left<n\\mathbf k|e^{-i\\mathbf{q\\cdot r}}|n'\\mathbf{k+q}\\right>|^2

    For ``gamma == 0``, this returns equation 17 of Ref 2, which is the
    polarizibility for general complex frequencies.

    For ``gamma > 0``, we return the Mermin-corrected Relaxation time approximation
    (Eqn 4.9 of Ref 1), which calls the gamma==0 case.


    Examples
    --------

    Plot the real and imaginary part of :math:`\\chi^0`, normalized to density of states. Replicates Fig. 1 of Ref. [3].

    .. plot::

        >>> from graphenemodeling.graphene import monolayer as mlg
        >>> import matplotlib.pyplot as plt
        >>> from matplotlib import cm
        >>> from scipy.constants import elementary_charge as eV
        >>> from scipy.constants import hbar
        >>> eF = 0.4*eV
        >>> gamma = 0
        >>> DOS = mlg.DensityOfStates(eF,model='LowEnergy')
        >>> kF = mlg.FermiWavenumber(eF,model='LowEnergy')
        >>> omega = np.linspace(0.01,2.5,num=250) / (hbar/eF)
        >>> q = np.linspace(0.01,2.5,num=250) * kF
        >>> pol = mlg.Polarizibility(q, omega[:,np.newaxis],gamma,eF)
        >>> fig, (re_ax, im_ax)  = plt.subplots(1,2,figsize=(11,4))
        >>> re_img = re_ax.imshow(  np.real(pol)/DOS,
        ...                         extent=(q[0]/kF,q[-1]/kF,hbar*omega[0]/eF,hbar*omega[-1]/eF),
        ...                         vmin=-2,vmax=1, cmap=cm.gray,
        ...                         origin = 'lower', aspect='auto'
        ...                       )
        <...
        >>> re_cb = fig.colorbar(re_img, ax = re_ax)
        >>> im_img = im_ax.imshow(  np.imag(pol)/DOS,
        ...                         extent=(q[0]/kF,q[-1]/kF,hbar*omega[0]/eF,hbar*omega[-1]/eF),
        ...                         vmin=-1,vmax=0, cmap=cm.gray,
        ...                         origin = 'lower', aspect='auto'
        ...                       )
        <...
        >>> im_cb = fig.colorbar(im_img, ax = im_ax)
        >>> re_ax.set_ylabel('$\\hbar\\omega$/$\\epsilon_F$',fontsize=14)
        Text...
        >>> re_ax.set_xlabel('$q/k_F$',fontsize=14)
        Text...
        >>> im_ax.set_ylabel('$\\hbar\\omega$/$\\epsilon_F$',fontsize=14)
        Text...
        >>> im_ax.set_xlabel('$q/k_F$',fontsize=14)
        Text...
        >>> plt.show()

    References
    ----------

    [1] Christensen Thesis 2017

    [2] Sernelius, B. Retarded interactions in graphene systems.

    [3] Wunsch, B., Stauber, T., Sols, F., and Guinea, F. (2006).
    Dynamical polarization of graphene at finite doping. New J. Phys. 8, 318–318.
    https://doi.org/10.1088%2F1367-2630%2F8%2F12%2F318.


    [4] Hwang, E.H., and Das Sarma, S. (2007).
    Dielectric function, screening, and plasmons in two-dimensional graphene.
    Phys. Rev. B 75, 205418.
    https://link.aps.org/doi/10.1103/PhysRevB.75.205418.


    '''

    if gamma==0 and T==0:
        '''
        Equation 17 of Ref 2
        '''

        prefactor = -DensityOfStates(eFermi,model='LowEnergy')

        kF = FermiWavenumber(eFermi, model='LowEnergy')

        x = q / (2*kF)
        zbar = sc.hbar*omega / (2*eFermi)

        f = lambda x,zbar: (np.arcsin( (1-zbar)/x) + np.arcsin( (1+zbar)/x )
                        - ((zbar-1)/x)*np.sqrt(1 - ((zbar-1)/x)**2 )
                        + ((zbar+1)/x)*np.sqrt(1 - ((zbar+1)/x)**2 ) )


        dd = 1 + (x**2 / (4*np.sqrt(x**2 - (zbar+1e-9*1j)**2 ))) * (sc.pi - f(x,zbar+1e-9*1j))

        return prefactor * dd

    elif gamma !=0:
        # Mermin-corrected Relaxation Time Approximation (Eqn 4.9 of Ref 1)
        pol_complex_arg = Polarizibility(q,omega+1j*gamma,0,eFermi,T=0)
        pol_0 = Polarizibility(q,0,0,eFermi,T=0)

        numerator = (1 + 1j*gamma/omega) * pol_complex_arg
        denominator = 1 + ( 1j*gamma/omega * pol_complex_arg / pol_0 )
        return numerator / denominator

def dPolarizibility(q,omega,gamma,eFermi,T,dvar,diff=1e-7):
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
        P = lambda w: np.real(Polarizibility(q,w,gamma,eFermi,T))
        wa, wb = omega*(1-diff), omega*(1+diff)
        return (P(wb)-P(wa))/(2*omega*diff)

    elif dvar == 'q':
        P = lambda qv: np.real(Polarizibility(qv,omega,gamma,eFermi,T))
        qa,qb = q*(1-diff), q*(1+diff)
        return (P(qb)-P(qa))/(2*q*diff)

######################
# Optical Properties #
######################

def OpticalConductivity(q,omega,gamma,eFermi,T,model=None):
    '''The diagonal conductivity of graphene :math:`\\sigma_{xx}`.

    Parameters
    ----------

    q:          array-like
                wavenumbers at which to evaluate the nonlocal conductivity.
                Choose q=0 for the LOCAL conductivity.

    omega:      array-like
                angular frequency

    eFermi:     scalar
                the Fermi energy (J)

    gamma:      scalar
                scattering rate

    T:          scalar
                Temperature

    model:      string
                Typically 'None', but for a specific model, specify it.


    Returns
    -------

    array-like
    conductivity at every value of omega

    Examples
    --------
    
    Plot the optical conductivity normalized to intrinsic conductivity :math:`\\sigma_0`.
    Replicates Fig. 4.4 in Ref. [1].

    .. plot ::

        >>> from graphenemodeling.graphene import monolayer as mlg
        >>> from scipy.constants import elementary_charge, hbar
        >>> from graphenemodeling.graphene._constants import sigma_0
        >>> import matplotlib.pyplot as plt
        >>> eF = 0.4 * elementary_charge
        >>> g = 0.012 * elementary_charge / hbar
        >>> w = np.linspace(0.01,3,num=150) / (hbar/eF)
        >>> s_0K_0g = mlg.OpticalConductivity(q=0,omega=w,gamma=0,eFermi=eF,T=0)
        >>> s_0K_12g = mlg.OpticalConductivity(q=0,omega=w,gamma=g,eFermi=eF,T=0.01)
        >>> s_300K_12g = mlg.OpticalConductivity(q=0,omega=w,gamma=g,eFermi=eF,T=300)
        >>> fig, (re_ax, im_ax) = plt.subplots(1,2,figsize=(11,4))
        >>> s_Re = np.real(s_0K_0g)
        >>> s_Im = np.imag(s_0K_0g)
        >>> re_ax.plot(w*hbar/eF,s_Re/sigma_0,'r',label='T=0, $\\hbar\\gamma$=0 meV')
        <...
        >>> im_ax.plot(w*hbar/eF,s_Im/sigma_0,'r',label='T=0, $\\hbar\\gamma$=0 meV')
        <...
        >>> s_Re = np.real(s_0K_12g)
        >>> s_Im = np.imag(s_0K_12g)
        >>> re_ax.plot(w*hbar/eF,s_Re/sigma_0,color='royalblue',label='T=0, $\\hbar\\gamma$=12 meV')
        <...
        >>> im_ax.plot(w*hbar/eF,s_Im/sigma_0,color='royalblue',label='T=0, $\\hbar\\gamma$=12 meV')
        <...
        >>> s_Re = np.real(s_300K_12g)
        >>> s_Im = np.imag(s_300K_12g)
        >>> re_ax.plot(w*hbar/eF,s_Re/sigma_0,color='green',label='T=300 K, $\\hbar\\gamma$=12 meV')
        <...
        >>> im_ax.plot(w*hbar/eF,s_Im/sigma_0,color='green',label='T=300 K, $\\hbar\\gamma$=12 meV')
        <...
        >>> re_ax.set_ylabel('Re[$\\sigma$]/$\\sigma_0$')
        >>> re_ax.set_xlabel('$\\hbar\\omega/E_F$')
        >>> re_ax.plot(w*hbar/eF,np.ones_like(w),'--',color='gray')
        >>> re_ax.set_ylim(0,2)
        >>> im_ax.set_ylabel('Im[$\\sigma$]/$\\sigma_0$')
        >>> im_ax.set_xlabel('$\\hbar\\omega/E_F$')
        >>> im_ax.plot(w*hbar/eF,np.zeros_like(w),'--',color='gray')
        >>> im_ax.set_ylim(-2,3)
        >>> plt.legend()
        >>> plt.show()


    Notes
    -----
    The optical conductivity :math:`\\overleftrightarrow{\\sigma}(q,\\omega)` 
    relates the surface current :math:`\\mathbf K(\\omega)` 
    to an applied electric field :math:`\\mathbf E(\\omega)`

    .. math::

        \\mathbf K(\\omega)=\\int \\overleftrightarrow\\sigma(q,\\omega)\\mathbf E(\\omega) dq
    

    Here, :math:`\\omega` refers to the frequency and :math:`q` refers to the scattering wavevector. 
    In many cases, :math:`\\overleftrightarrow{\\sigma}` is isotropic, 
    so the above equation can be reduced to a scalar equation

    .. math::

        K(\\omega)=\\int \\sigma(q,\\omega)E(\\omega)dq


    The most general expression for the conductivity is given by the Kubo formula
    
    .. math::
        \\sigma(q,\\omega)=\\frac{ie^2\\omega}{q^2}\\chi^0(q,\\omega).
    
    If ``q`` is nonzero, this form is used.
    However, it is common to use simpler limiting cases of this expression.
    The local conductivity (``q=0``) is the one which is most familiar 
    and it relates the surface current to the electric field linearly
    
    .. math::

        \\mathbf K(\\omega)=\\sigma(\\omega)\\mathbf E

    It can be found from the nonlocal conductivity by taking the limit $\\lim_{q\\to 0}\\sigma(q,\\omega)=\\sigma(\\omega)$.

    References
    ----------

    [1] Christensen, T. (2017).
    From Classical to Quantum Plasmonics in Three and Two Dimensions (Cham: Springer International Publishing).
    http://link.springer.com/10.1007/978-3-319-48562-1


    '''

    # Local case
    if np.all(q) == 0:

        if T!=0:
            intra_pre = 4 * _c.sigma_0 * (2*1j*sc.k*T) / (sc.pi*sc.hbar)
            inter_pre = _c.sigma_0

            ### Intraband Contribution ###

            # Using np.logaddexp() avoids the large cosh in ln( cosh(1/T) )
            x = eFermi / (2*sc.k*T)
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
            intra = lambda w: 1j*_c.sigma_0 * 4*eFermi / (sc.pi*sc.hbar* (w + 1j*gamma))

            inter = lambda w: _c.sigma_0 * ( np.heaviside(sc.hbar*w - 2*eFermi,0.5) + 
                                            (1j/sc.pi) * np.log(np.abs((2*eFermi-sc.hbar*w)/(2*eFermi+sc.hbar*w))))

            conductivity = intra(omega) + inter(omega)

    # Nonlocal case
    else:
        if T==0:
            conductivity = 1j*sc.elementary_charge**2 * (omega / q**2) * Polarizibility(q,omega,gamma,eFermi,T)

        else:
            pass

    return conductivity

def OpticalConductivityMatrix(q,omega,gamma, eFermi,T,mu0,mu0T):
    '''
    Returns the conductivity matrix of monolayer graphene.

    Parameters
    ----------

    omega:

    Returns
    ----------

    sigma_matrix:   2x2 numpy array, conductivity of monolayer graphene

    '''


    conductivity_matrix = np.array([[OpticalConductivity(q,omega,gamma,eFermi,T),0],
                                    [0,OpticalConductivity(q,omega,gamma,eFermi,T)]])

    return conductivity_matrix

def Permittivity(omega,eFermi,T, gamma=None,epsR=None,model=None):
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
        x4 = _c.vF
        x5 = Mobility(T,mu0,mu0T) # mobility at the new temperature
        x6 = epsR

        x7 = _c.thickness # 3.4 angstroms by default
        x8 = sc.epsilon_0

        prefactor = x1**2*x3 / ( sc.pi * x2**2)
        denominator = omega**2 + ( x1*x4**2 / (x5*x3) )**2

        term1 = - prefactor*(omega*x8*x7)**(-1) * omega / denominator
        term2 = 1j*prefactor * (x1*x4**2 / (omega*x5*x3*x8*x7)) / denominator

        eps = x6 + term1 + term2

        return eps

    else:
        eps = 1 + 1j*OpticalConductivity(0,omega,gamma,eFermi,T)/(omega*sc.epsilon_0)

    return eps

def FresnelReflection(kpar,omega,gamma,eFermi,T,eps1,eps2,polarization):
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

    sigma = OpticalConductivity(kpar,omega,gamma,eFermi,T)

    if polarization=='p' or polarization=='TM':
        numerator   = eps2*kperp1 - eps1*kperp2 + ( sigma*kperp1*kperp2 / (sc.epsilon_0*omega) )
        denominator = eps2*kperp1 + eps1*kperp2 + ( sigma*kperp1*kperp2 / (sc.epsilon_0*omega) )

    if polarization=='s' or polarization=='TE':
        numerator = kperp1 - kperp2 - sc.mu_0*omega*sigma
        denominator = kperp1 + kperp2 + sc.mu_0*omega*sigma

    return numerator / denominator

def LocalDensityOfStates(d,omega,gamma,eFermi,T):
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
        Im_rp      = lambda kpar: np.imag( FresnelReflection(kpar,omega,gamma,eFermi,T,1,1,'p') )

        integrand = lambda kpar: (kpar**2/k0w**3)*Im_rp(kpar)*np.exp(-2*kpar*d0)

        a,b = 1e-3, np.abs(_c.K) # Increasing this bound does not lead to more convergence

        k_plasmon=InversePlasmonDispersion(omega,gamma,eFermi,1,1,T,model='nonlocal')

        integral[i] = integrate.quad(integrand,a,b,
                                    points=(k_plasmon),limit=100)[0]

    return ldos0 * integral      

def OpticalResponseBound(omega,gamma,eFermi,T,d=None,method='svd',restype='material'):
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
        sigma = OpticalConductivity(0,omega,gamma,eFermi,T)

        bound = Z_0 * np.abs(sigma)**2 / np.real(sigma)

    elif method == 'svd':
        bound = np.empty_like(omega)

        for i, w in np.ndenumerate(omega):
            s = OpticalConductivityMatrix(0,w,gamma,eFermi,T)
            s_dag = np.conj(np.transpose(s))
            s_re_inv = np.linalg.inv(np.real(s))
            prod = np.dot(s_dag,np.dot(s_re_inv,s))
            two_norm= np.linalg.svd(prod)[1][0]
            bound[i] = Z_0*two_norm



    return prop[restype] * bound

#####################
# Phonon Properties #
#####################

def PhononSelfEnergy(self):
    pass

##############
# Plasmonics #
##############

def PlasmonDispersion(q,gamma,eFermi,eps1,eps2,T,model):
    '''
    Returns the nonretarded plasmon dispersion relations E(q) for a surface
    plasmon in an infinite sheet of graphene sandwiched between two
    dielectrics eps1 and eps2.

    All values returned are assumed to be at zero temperature with no loss (gamma).

    Parameters
    ----------

    q:  array-like
        Wavenumber of the plasmon

    eps1,eps2:  scalar
        the relative permittivity of upper and lower half-space respectively.

    model:  string
        'intra' for intraband dispersion,
        'local' uses the intraband + interband constributions to the conductivity, and
        'nonlocal' for fully nonlocal conductivity.

    Returns
    -------

    omega:  array-like
        Frequency of the plasmon with wavenumber q

    Notes
    -----

    ``model=='intra'`` uses

    .. math::

        \\omega = \\frac{1}{\\hbar}\\sqrt{\\frac{e^2\\epsilon_F}{2\\pi\\epsilon_0\\bar\\epsilon}q}

    ``model=='local'`` uses

    .. math::

        1-\\frac{i\\text{Im}[\\sigma(q,\\omega)]}{2i\\epsilon_0\\bar\\epsilon\\omega}=0

    and finally, ``model=='nonlocal'`` uses

    .. math::

        1 - \\frac{\\sigma(q,\\omega)q}{2i\\epsilon_0\\bar\\epsilon\\omega} = 0


    Examples
    --------

    Plot the three expressions for the dispersion relation. Replicates Fig. 5.2 in Ref. [1].

    .. plot::

        >>> from graphenemodeling.graphene import monolayer as mlg
        >>> from scipy.constants import elementary_charge, hbar
        >>> eV = elementary_charge
        >>> gamma=0.012 * eV / hbar
        >>> eF = 0.4*eV
        >>> kF = mlg.FermiWavenumber(eF,model='LowEnergy')
        >>> q = np.linspace(1e-3,3,num=200) * kF
        >>> disp_intra = mlg.PlasmonDispersion(q,gamma,eF,eps1=1,eps2=1,T=0,model='intra')
        >>> disp_local = mlg.PlasmonDispersion(q,gamma,eF,eps1=1,eps2=1,T=0,model='local')
        >>> disp_nonlocal = mlg.PlasmonDispersion(q,gamma,eF,eps1=1,eps2=1,T=0,model='nonlocal')
        >>> fig, ax = plt.subplots(figsize=(6,6))
        >>> ax.plot(q/kF,hbar*disp_intra/eF,'--',label='Intraband')
        <...
        >>> ax.plot(q/kF,hbar*disp_local/eF,'r-.',label='Local')
        <...
        >>> ax.plot(q[:190]/kF,hbar*disp_nonlocal[:190]/eF,'g',label='Nonlocal')
        <...
        >>> ax.plot(q/kF,q/kF,color='gray',linestyle='--')
        <...
        >>> ax.set_xlabel('$q/k_F$')
        >>> ax.set_ylabel('$\\hbar\\omega/E_F$')
        >>> ax.set_xlim(0,3)
        >>> ax.set_ylim(0,3)
        >>> plt.legend()
        >>> plt.show()

    Plot dispersion relation with a lower half-space permittivity of :math:`\\epsilon=4`` (an approximation for hexagonal boron nitride).
    Replicates Fig. 1d in Ref. [2].

    .. plot::

        >>> from graphenemodeling.graphene import monolayer as mlg
        >>> from scipy.constants import elementary_charge, hbar
        >>> eV = elementary_charge
        >>> gamma=0.012 * eV / hbar
        >>> eF = 0.4*eV
        >>> kF = mlg.FermiWavenumber(eF,model='LowEnergy')
        >>> q = np.linspace(1e-3,2,num=200) * kF
        >>> vac_intra = mlg.PlasmonDispersion(q,gamma,eF,eps1=1,eps2=1,T=0,model='intra')
        >>> hbn_intra = mlg.PlasmonDispersion(q,gamma,eF,eps1=1,eps2=4,T=0,model='intra')
        >>> vac_nonlocal = mlg.PlasmonDispersion(q,gamma,eF,eps1=1,eps2=1,T=0,model='nonlocal')
        >>> hbn_nonlocal = mlg.PlasmonDispersion(q,gamma,eF,eps1=1,eps2=4,T=0,model='nonlocal')
        >>> fig, ax = plt.subplots(figsize=(6,6))
        >>> ax.plot(q/kF,hbar*vac_intra/eF,'k-.',label='Vacuum Intraband')
        <...
        >>> ax.plot(q/kF,hbar*hbn_intra/eF,linestyle='dotted',color='purple',label='hBN Intraband')
        <...
        >>> ax.plot(q/kF,hbar*vac_nonlocal/eF,'k-',label='Vacuum Nonlocal')
        <...
        >>> ax.plot(q/kF,hbar*hbn_nonlocal/eF,'--',color='purple',label='hBN Nonlocal')
        <...
        >>> ax.plot(q/kF,q/kF,color='gray',linestyle='--')
        <...
        >>> ax.set_xlabel('$q/k_F$')
        >>> ax.set_ylabel('$\\hbar\\omega/E_F$')
        >>> ax.set_xlim(0,2)
        >>> ax.set_ylim(0,2)
        >>> plt.legend()
        >>> plt.show()

    References
    ----------

    [1] Christensen, T. (2017).
    From Classical to Quantum Plasmonics in Three and Two Dimensions (Cham: Springer International Publishing).
    http://link.springer.com/10.1007/978-3-319-48562-1.

    [2] Jablan, M., Buljan, H., and Soljačić, M. (2009).
    Plasmonics in graphene at infrared frequencies. Phys. Rev. B 80, 245435.
    https://link.aps.org/doi/10.1103/PhysRevB.80.245435.


    '''

    epsavg = (eps1+eps2)/2

    # Analytical expression in intraband approximation
    if model=='intra':
        radical = q * sc.elementary_charge**2 * eFermi / (2*sc.pi * sc.epsilon_0 * epsavg)
        return (1/sc.hbar) * np.sqrt(radical)

    if model=='local':
        omega = np.empty_like(q)

        sigma = lambda w: OpticalConductivity(0,w,gamma,eFermi,T=0)

        for i,q0 in np.ndenumerate(q):
            root_eqn = lambda w: 1 - np.imag(sigma(w))*q0 / (2*sc.epsilon_0*epsavg*w)

            a, b = PlasmonDispersion(q0,gamma,eFermi,eps1,eps2,T,model='intra'), 1e-8
            omega[i] = optimize.brentq(root_eqn,a,b)

        return omega

    if model=='nonlocal':
        omega = np.empty_like(q)

        kF = FermiWavenumber(eFermi,model='LowEnergy')

        for i, q0 in np.ndenumerate(q):
            root_eqn = lambda w: PlasmonDispersionRoot(q0,w,gamma,eFermi,eps1,eps2,T=0)
            
            # Frequency is definitely below 1,1 intraband dispersion
            b = PlasmonDispersion(q0,gamma,eFermi,1,1,T,model='intra')

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

def PlasmonDispersionRes(q,gamma,eFermi,eps1,eps2,T,exp_res=1):
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
        w1 = q0*_c.vF
        w2 = PlasmonDispersion(q0,gamma,eFermi,eps1,eps2,T,model='intra')
        w0 = (w2+w1)/2
        # omega=np.linspace(q0*_c.vF,2*q0*_c.vF,num=300)
        omega=np.linspace(w1,w2,num=300)

        y = np.imag(FresnelReflection(q0,omega,gamma,eFermi,T,eps1,eps2,'TM'))

        p0=[w0,0.1*w0,10]

        fit = optimize.leastsq(ErrorFunc,p0,
                                    args=(omega,y))

        pFit[i,:] = fit[0]



    return np.abs(pFit)

def InversePlasmonDispersion(omega,gamma,eFermi,eps1,eps2,T,model):
    '''
    Returns the wavenumber of a plasmon given the frequency.

    Useful when determining the wavelength of a plasmon excited by light.
    '''

    kF = FermiWavenumber(eFermi,model='LowEnergy')

    cutoff = 4*kF
    q = np.empty_like(omega)

    for i, omega in np.ndenumerate(omega):
        root_eqn = lambda q: np.abs( omega - PlasmonDispersion(q,gamma,eFermi,eps1,eps2,T,model) )

        reps = 1
        while reps < 5:
            q[i] = optimize.minimize_scalar(root_eqn,bounds=(1e-6,reps*cutoff),method='bounded').x

            if q[i] >= cutoff:
                reps=reps+1
            else:
                reps=5

    return q

def PlasmonDispersionRoot(q,omega,gamma,eFermi, eps1,eps2 ,T):
    '''
    The equation used for numerically solving the plasmon dispersion in the nonretarded regime.
    '''

    epsavg = (eps1+eps2)/2

    return 1 - np.imag(OpticalConductivity(q,omega,gamma,eFermi,T))*q / (2*sc.epsilon_0*epsavg*omega)

def PlasmonDispersionLoss(omega,gamma,eFermi,eps1,eps2,T,model):
    '''
    The loss of the plasmon wavenumber q=q1+iq2. Returns q2. Equation 15 of Ref [1]
    with tau = infinity (or gamma = 0). Assumes q2<<q1

    References
    ----------

    [1] Jablan et al. 2009 (note tau = 1/ gamma)
    '''
    
    if model == 'nonlocal':
        q1 = InversePlasmonDispersion(omega,gamma,eFermi,eps1,eps2,T,model)

        pol = Polarizibility(q1,omega,gamma,eFermi,T=0)
        pol0= Polarizibility(q1,1e-9,gamma,eFermi,T=0)

        dpolq = dPolarizibility(q1,omega,gamma,eFermi,T,dvar='q')
        dpolw = dPolarizibility(q1,omega,gamma,eFermi,T,dvar='omega')

        numerator = np.imag(-pol) + gamma * (-1)*dpolw + (gamma/omega) * np.real(-pol * (1- (pol/pol0)))
        denominator = (1/q1) * np.real(-pol) - (-1)*dpolq

        q2 = numerator / denominator

    return q2

def dPlasmonDispersion(q,gamma,eFermi,eps1,eps2,T,model,dvar=None,diff=1e-7):
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
            w = lambda eF: PlasmonDispersion(q0,gamma,eF,eps1,eps2,T,model)
            e1, e2 = eFermi*(1-diff), eFermi*(1+diff)
            result[i] = (w(e2)-w(e1))/(2*eFermi*diff)

        return result

    if dvar=='q':
        pass

    pass

def DipoleDecayRate(z,omega,gamma,eFermi,T,eps1,eps2):
    '''
    Decay rate of a dipole emitter placed near graphene.
    Right now, the dipole only points perpendicular.

    This should be moved to the Emitter.Dipole object

    Eqn 5 in the SM of Ref 1.

    References
    ----------

    [1] Koppens et al. 2011
        URL: https://doi.org/10.1021/nl201771h
    '''
    warnings.warn('Monolayer.DipoleDecayRate: Ignoring dipole components in xy-plane')

    dipole = Emitter.Dipole()

    d = np.array([0,0,1])

    integral = np.empty_like(omega)

    for i, w in np.ndenumerate(omega):

        kperp   = lambda kpar: np.sqrt(eps1*(w/sc.speed_of_light)**2 - kpar**2)
        rp      = lambda kpar: FresnelReflection(kpar,w,gamma,eFermi,T,eps1,eps2,'p')
        perpterm = lambda kpar: 2*np.abs(d[2])**2 * kpar**2 * rp(kpar)

        integrand = lambda kpar: np.real( (kpar - 1e-9*1j) * np.real( perpterm(kpar-1e-9*1j) * np.exp(2*1j*kperp(kpar-1e-9*1j)*z) / kperp(kpar-1e-9*1j) ) )

        b = np.abs(_c.K) # Increasing this bound does not lead to more convergence
        kpar_pol = np.sqrt(eps1) * (w/sc.speed_of_light)
        k_plasmon=InversePlasmonDispersion(w,gamma,eFermi,eps1,eps2,T,model='nonlocal')

        integral[i] = integrate.quad(integrand,1e-3,b,
                                    points=(kpar_pol,k_plasmon),limit=100)[0]

    return dipole.DecayRate(omega,d) + (1/sc.hbar) * integral
