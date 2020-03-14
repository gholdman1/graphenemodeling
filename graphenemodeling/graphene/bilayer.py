"""
======================================================
Bilayer (:mod:`graphenemodeling.graphene.bilayer`)
======================================================

Functions
=========

Band Structure
--------------

.. toctree::
    :maxdepth: 1

    graphene.bilayer.Hamiltonian

"""

import numpy as np
import scipy.constants as sc

from graphenemodeling.graphene.base import BaseGraphene
eVtoJ = sc.elementary_charge
e0 = sc.epsilon_0

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

######################
### Band Structure ###
######################

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

def CarrierDispersion(self,k,u,band,approx='Common'):
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

def DensityOfStates(self, e, u):
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

def FermiLevel(self,n,u):
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

    vplus:  Scalar; the potential of the BLG. Related to fermi level by vplus = eF / elementary_charge

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
