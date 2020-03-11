import os
import numpy as np
import warnings

from scipy import integrate
from scipy import optimize

import scipy.constants as sc

def constant_to_callable(eps):
	'''
	Takes in a constant eps (int, float, or complex) and returns a function
	eps(kpar,omega).

	Useful for writing code that assumes callables for generality.
	'''

	epsval=eps
	del eps
	eps = lambda kp, w: epsval

	return eps

def ScalarGreenFunctionVacuum(r1,r2,omega):
	'''
	The scalar Green Function in vacuum

	See Ref 1

	Parameters
	----------

	r1,r2:		length 3 arrays, the two positions at which to evaluate the Green function.

	omega:		scalar, frequency at which to evaluate the Green function.

	References
	----------

	[1] Sarabandi, K. n.d. “Dyadic Green’s Functions.”
		https://www.eecs.umich.edu/courses/eecs730/lect/DyadicGF_W09_port.pdf.

	'''

	k0 = omega/sc.speed_of_light
	dr = r1-r2

	if np.dot(dr,dr) == 0:
		return 1j*k0 / (4*np.pi)

	absdr = np.sqrt(np.dot(dr,dr))

	num = np.exp(1j*k0*absdr)
	den = 4*np.pi*absdr

	G = num/den 

	return G

def DyadicGreenFunctionVacuum(r1,r2,omega):
	'''
	The dyadic Green function in vacuum.	

	See Ref 1

	Parameters
	----------

	r1:			length 3 array, one of the positions at which to evaluate the Green function.

	r2:			length 3 array, see r1

	omega:		array-like, frequency at which to evaluate the Green function.

	References
	----------

	[1] Sarabandi, K. n.d. “Dyadic Green’s Functions.”
		https://www.eecs.umich.edu/courses/eecs730/lect/DyadicGF_W09_port.pdf.

	'''

	dr = r1-r2
	absdr=np.sqrt(np.dot(dr,dr))
	RR = np.outer(dr/absdr,dr/absdr)
	I = np.eye(3)

	k0 = omega/sc.speed_of_light

	term1 = (3 / (k0*absdr)**2) - (3*1j / (k0*absdr) ) - 1
	term2 = 1 + (1j / (k0*absdr)) - (1 / (k0*absdr)**2)

	G = ( term1 * RR + term2 * I ) * ScalarGreenFunctionVacuum(r1,r2,omega)

	return G

def ScalarGreenFunctionHalfSpace(r1,r2,omega,eps1,eps2):
	'''

	Parameters
	---------

	eps1:	permittivity in z>0

	eps2:	permittivity in z<0

	omega:	

	References
	-----------

	[1] Barcellona, Pablo, Robert Bennett, and Stefan Yoshi Buhmann. 2018.
	“Manipulating the Coulomb Interaction: A Green\textquotesingles Function Perspective.”
	Journal of Physics Communications 2 (3): 035027. https://doi.org/10.1088/2399-6528/aaa70a.

	'''
	warnings.warn('ScalarGreenFunctionHalfSpace: Not tested.\nMay only be valid for real, constant eps1 and eps2.\nAlso does not include conductive interface.')
	
	r2im = np.copy(r2)
	r2im[2] = -r2[2] # Image charge of source

	dr = r1-r2
	drim = r1-r2im

	absdr = np.sqrt(np.dot(dr,dr))
	absdrim=np.sqrt(np.dot(drim,drim))

	k = 1 / (4*np.pi)

	# Source is in eps1
	if r2[2] > 0:
		if r1[2]>0:
			G = (eps1**(-1)) * ( (k/absdr) + (k/absdrim)*(eps1-eps2)/(eps1+eps2) )
		else:
			G = ( 2 / (eps1+eps2) ) * (k*absdr)

	# Source is in eps2, swap eps1<->eps2 and condition on r1
	if r2[2] < 0:
		if r1[2]<0:
			G = (eps2**(-1)) * ( (k/absdr) + (k/absdrim)*(eps2-eps1)/(eps2+eps1) )
		else:
			G = ( 2 / (eps2+eps1) ) * (k*absdr)


	return G

def DyadicGreenFunctionHalfSpace(r1,r2,omega,eps1,eps2):
	pass


class GreenFunction:
	'''
	A class for defining Green Functions in various spaces

	'''


	def __init__(self,field,space,eps=None):
		self.field=field
		self.space=space

		self.eps=eps # permittivity of lower half space


	def Scalar(self,r1,r2,omega):
		'''
		The scalar Greens function in the space
		'''

		if self.space=='vacuum':
			return ScalarGreenFunctionVacuum(r1,r2,omega)

		if self.space=='halfspace':
			return ScalarGreenFunctionHalfSpace(r1,r2,omega,eps)
		else:
			return None

	def Dyad(self,r1,r2,omega):
		'''


		Parameters
		----------

		omega:		array-like, frequency at which to evaluate (rad/s)
		'''
		if self.space=='vacuum':
			if self.field=='electric':
				return DyadicGreenFunctionVacuum(r1,r2,omega)
			if self.field=='magnetic':
				return DyadicGreenFunctionVacuum(r1,r2,omega)
		
		if self.space=='halfspace':
			if self.field=='electric':
				return DyadicGreenFunctionHalfSpace(r1,r2,omega,eps)
		else:
			return None

	def Trace(self,r1,r2,omega):
		'''
		Trace of self.Dyad. May be necessary if Dyad is singular at r1=r2.

		Parameters
		----------

		omega:	array-like, frequency at which to evaluate (rad/s)
		'''

		if self.space=='vacuum':
			return 2*ScalarGreenFunctionVacuum(r1,r2,omega)

def FresnelReflection(kpar,omega,eps1,eps2,pol,sigma):
	'''
	Returns the Fresnel reflection coefficient from a boundary between dielectrics.
	Uses equation 2.49 of Ref [1]

	Light is incident from medium 1 (eps1) onto medium 2 (eps2)

	Assumes non-magnetic materials.

	Parameters
	----------

	kpar:		array-like, the in-plane momentum

	omega:		array-like, the frequency on incident light

	eps1:		scalar or callable, the relative permittivity of the 
				half space from which light is incident.

	eps2:		scalar or callable, the relative permittivity of the 
				half space to which light is transmitted.

	pol:		polarization, 'p' ('TM') or 's' ('TE')

	sigma:		float or callable, the conductivity of a 2D material
				between the medium1 and medium 2.s

	References
	----------

	[1] Novotny, Lukas, and Bert Hecht. 2006. Principles of Nano-Optics.
		Cambridge University Press.
		http://www.fulviofrisone.com/attachments/article/406/Principles%20of%20Nano-Optics.pdf.

	[2] Christensen, Thomas. 2017. From Classical to Quantum Plasmonics in Three and Two Dimensions.
		Springer Theses. Cham: Springer International Publishing.
		https://doi.org/10.1007/978-3-319-48562-1.
	'''

	if not callable(eps1):
		eps1float=eps1
		del eps1
		eps1 = lambda kp, w: eps1float
	if not callable(eps2):
		eps2float=eps2
		del eps2
		eps2 = lambda kp, w: eps2float
	if not callable(sigma):
		sigfloat = sigma
		del sigma
		sigma = lambda kp,w: sigfloat

	k0 = (omega/sc.speed_of_light)

	# Get perpendicular momenta
	kperp1 = np.sqrt(eps1(kpar,omega)*k0**2 - kpar**2 + 1e-9*1j)
	kperp2 = np.sqrt(eps2(kpar,omega)*k0**2 - kpar**2 + 1e-9*1j)

	if pol =='p' or pol=='TM':

		numerator = (eps2(kpar,omega)*kperp1
					-eps1(kpar,omega)*kperp2
					+ sigma(kpar,omega)*kperp1*kperp2/(sc.epsilon_0*omega) )
		denominator= (eps2(kpar,omega)*kperp1
					+eps1(kpar,omega)*kperp2
					+ sigma(kpar,omega)*kperp1*kperp2/(sc.epsilon_0*omega) )

		return numerator / denominator

	if pol=='s' or pol=='TE':
		numerator = kperp1-kperp2 - sc.mu_0*omega*sigma(kpar,omega)
		denominator=kperp1+kperp2 + sc.mu_0*omega*sigma(kpar,omega)

		return numerator/denominator

def FresnelTransmission(kpar,omega,eps1,eps2,pol,sigma=0):
	'''
	Returns the Fresnel reflection coefficient from a boundary between dielectrics.
	Uses Eqn 2.50 of Ref [1]

	Light is incident from medium 1 (eps1) on medium 2 (eps2)

	Currently only p polarization.

	Assumes non-magnetic materials.

	Parameters
	----------

	kpar:		array-like, the in-plane momentum

	omega:		array-like, the frequency on incident light

	eps1:		scalar or callable, the relative permittivity of the 
				half space from which light is incident.

	eps2:		scalar or callable, the relative permittivity of the 
				half space to which light is transmitted.

	pol:		polarization, currently only 'p' ('TM')

	References
	----------

	[1] Novotny, Lukas, and Bert Hecht. 2006. Principles of Nano-Optics.
		Cambridge University Press.
		http://www.fulviofrisone.com/attachments/article/406/Principles%20of%20Nano-Optics.pdf.

	[2] Christensen, Thomas. 2017. From Classical to Quantum Plasmonics in Three and Two Dimensions.
		Springer Theses. Cham: Springer International Publishing.
		https://doi.org/10.1007/978-3-319-48562-1.

	'''

	# Define constant permittivities as callables
	# for generality.
	if not callable(eps1):
		eps1float=eps1
		del eps1
		eps1 = lambda kp, w: eps1float
	if not callable(eps2):
		eps2float=eps2
		del eps2
		eps2 = lambda kp, w: eps2float
	if not callable(sigma):
		sigfloat = sigma
		del sigma
		sigma = lambda kp,w: sigfloat

	k0 = (omega/sc.speed_of_light)
	# Get perpendicular momenta
	kperp1 = np.sqrt(eps1(kpar,omega)*k0**2 - kpar**2 + 1j*1e-9)
	kperp2 = np.sqrt(eps2(kpar,omega)*k0**2 - kpar**2 + 1j*1e-9)

	if pol =='p' or pol=='TM':

		numerator = 2 * eps2(kpar,omega)*kperp1
		denominator= ( eps2(kpar,omega)*kperp1 + 
					   eps1(kpar,omega)*kperp2 +
					   sigma(kpar,omega)*kperp1*kperp2 / (sc.epsilon_0*omega) )

		return (numerator / denominator) * np.sqrt(eps1(kpar,omega)/eps2(kpar,omega))

	if pol=='s' or pol=='TE':

		numerator=2*kperp1
		denominator=kperp1 + kperp2 + sc.mu_0*omega*sigma(kpar,omega)

		return numerator/denominator

def FresnelTransferMatrix(kpar,omega,eps1,eps2,pol,sigma):
	'''
	Transfer matrix relating the electric fields in eps1 and eps2 across an interface.

	E1 = FresnelTransferMatrix E2

	where Ei is a vector giving the electric field propagating in the +z and -z directions respectively.

	Equation 3.10 of Reference 1

	Parameters
	----------

	kpar:		scalar, the in-plane momentum

	omega:		scalar, the frequency

	eps1:		scalar or callable, the relative permittivity of the 
				half space from which light is incident.

	eps2:		scalar or callable, the relative permittivity of the 
				half space to which light is transmitted.

	pol:		polarization, 'p' ('TM') or 's' ('TE')

	sigma:		float or callable, the conductivity of a 2D material
				between the medium1 and medium 2.s				

	Returns
	----------

	FTM:		2 x 2 array, the Fresnel Transfer Matrix

	References
	----------

	[1] Sipe, J. E. 1987. “New Green-Function Formalism for Surface Optics.”
		Journal of the Optical Society of America B 4 (4): 481.
		https://doi.org/10.1364/JOSAB.4.000481.

	'''

	r12 = FresnelReflection(kpar,omega,eps1,eps2,pol,sigma)
	t12 = FresnelTransmission(kpar,omega,eps1,eps2,pol,sigma)

	prefactor = 1 / t12

	matrix = np.array([[1,r12],[r12,1]])

	FTM = prefactor*matrix

	return FTM

def PropagationMatrix(z,kpar,omega,eps):
	'''
	Propagates electric field a distance z in the positive z direction.

	Equation 3.15 of Ref 1

	References
	----------

	[1] Sipe, J. E. 1987. “New Green-Function Formalism for Surface Optics.”
		Journal of the Optical Society of America B 4 (4): 481.
		https://doi.org/10.1364/JOSAB.4.000481.

	'''

	k0 = omega / sc.speed_of_light

	# Define constant permittivities as callables
	# for generality.
	if not callable(eps):
		epsfloat=eps
		del eps
		eps = lambda kp, w: epsfloat

	kperp = np.sqrt(eps(kpar,omega)*k0**2 - kpar**2 + 1e-9*1j)

	matrix = np.array([[np.exp(1j*kperp*z),0],[0,np.exp(-1j*kperp*z)]])

	return matrix

def FresnelTransferMatrixStack(kpar,omega,eps,sigma,d,pol):
	'''

	Generalization of equations 3.17 in Ref1

	Parameters
	----------

	kpar:		in-plane momentum

	omega:		frequency

	eps:		array of shape (n+2) where n is the number of finite dielectrics.
				Each entry is given by [epsi] where epsi is a constant or callable
				giving the permittivity and di is the thickness of that layer.
				First and last thickness must be -1.

	sigma:		array of shape (n+1). Conductivity of interfaces
				sigma[0] is the conductivity of the interface between eps[0] and eps[1]

	d:			array of shape n+2 giving the thickness of dielectrics 0 through n+1.
				First and last entries must be -1

	References
	----------

	[1] Sipe, J. E. 1987. “New Green-Function Formalism for Surface Optics.”
		Journal of the Optical Society of America B 4 (4): 481.
		https://doi.org/10.1364/JOSAB.4.000481.
	'''

	# Error checking
	if d[0]!=-1 or d[-1]!=-1:
		raise ValueError('First and last entries of d must be -1')
	if np.shape(eps)[0]-1 != np.shape(sigma)[0]:
		raise ValueError('First axis of eps must have length exactly one larger than first axis of sigma')

	# Number of finite layers
	N = np.shape(eps)[0]-2

	# Initial transfer matrix
	M = np.eye(2)

	for i in range(N):
		j = N-i

		# Cross interface
		M = np.dot(FresnelTransferMatrix(kpar,omega,eps[j],eps[j+1],pol,sigma[j]),M)

		# Propagate to the next interface
		M = np.dot(PropagationMatrix(d[j],kpar,omega,eps[j]),M)


	M = np.dot(FresnelTransferMatrix(kpar,omega,eps[0],eps[1],pol,sigma[0]),M)

	return M

def FresnelCoefficientStack(kpar,omega,eps,sigma,d,pol):

	R = np.empty((np.size(kpar),np.size(omega)),dtype=complex)
	T = np.empty_like(R)

	for i in range(np.size(kpar)):
		for j in range(np.size(omega)):
			M = FresnelTransferMatrixStack(kpar[i],omega[j],eps,sigma,d,pol)
			T[i,j] = 1 / M[1,1]
			R[i,j] = M[0,1] / M[1,1]

	return (R,T)



def PlasmonDispersionRelation(kpar,eps1,eps2,sigma=0,
							start=lambda x: [1,1],
							valid_only=False,display=False):
	'''
	The plasmon dispersion relation, AKA the poles of
	FresnelReflection.
	Computed by finding the complex omega which is a 
	zero of the denominator

	Parameters
	----------

	kpar:		array-like, in-plane momentum (1/m)

	eps1:		number or callable, permittivity in z>0

	eps2:		number or callable, permittivity in z<0

	sigma:		number or callable, conductivity of the interface.

	start:		callable returning [omega1, omega2] where omega1 is a frequency
				near the desired branch of the dispersion relation.
				omega2 is an estimate of the loss.

	validonly:	Boolean, Choose "False" to return all omega which were found via a fit.
						Choose "True" to return tuple (kpar, omega) which return a true root

	display:	Boolean, whether to display the output of fmin
	'''

	kpar = np.atleast_1d(kpar)
	# Define constant permittivities as callables.
	if type(eps1)==float or type(eps1)==int:
		eps1float=eps1
		del eps1
		eps1 = lambda kp, w: eps1float
	if type(eps2)==float or type(eps2)==int:
		eps2float=eps2
		del eps2
		eps2 = lambda kp, w: eps2float
	if type(sigma)==float or type(sigma)==int:
		sigfloat = sigma
		del sigma
		sigma = lambda kp,w: sigfloat

	# Get perpendicular momenta
	kperp1 = lambda kp,w: np.sqrt( (kp*sc.speed_of_light/w)**2 - eps1(kp,w) - 1e-9*1j)
	kperp2 = lambda kp,w: np.sqrt( (kp*sc.speed_of_light/w)**2 - eps2(kp,w) - 1e-9*1j)

	term1 = lambda kp,w: eps1(kp,w) / kperp1(kp,w)
	term2 = lambda kp,w: eps2(kp,w) / kperp2(kp,w)
	term3 = lambda kp,w: sigma(kp,w) / (1j*sc.epsilon_0*sc.speed_of_light)

	# Initialize as [[-1,-1],...]
	# [-1,-1] is a signal of invalid solution
	omega1 = np.zeros((np.size(kpar),2)) -1

	startval = start(kpar[0])
	for i,kp in np.ndenumerate(kpar):

		# Define plasmon dispersion relation (pdr) of which to find zero
		pdr = lambda w:np.abs( term1(kp,w[0]-1j*w[1]) + term2(kp,w[0]-1j*w[1]) - term3(kp,w[0]-1j*w[1]) ) 

		result = optimize.fmin(pdr,startval,disp=display)

		# Should check that it's actually a zero
		
		if valid_only == True:
			if pdr(result)<1e-2:
				# If successful, record result and use as starting point for next iteration.
				omega1[i[0]] = result
				startval = omega1[i[0]]
			else:
				# Otherwise, signal it failed and start from suggested startval
				omega1[i[0]] = [-1,-1]
				startval = start(kpar[i[0]])
		else:
			omega1[i[0]] = result
			startval = omega1[i[0]]


	if valid_only == True:
		omegares1 = omega1[omega1[:,0]>0,:]
		kpar1 = kpar[omega1[:,0]>0]
		return (kpar1,omegares1)

	return (kpar, omega1)


def InversePlasmonDispersionRelation(omega,eps1,eps2,sigma=0,
									start=lambda x: [1,1],qmin=1e5,qmax=1e9,qnum=100):
	'''		
	Finds value of in-plane momentum that yields omega. Assumes monotonic dispersion.
	'''
	omega=np.atleast_1d(omega)

	qreturn = np.empty_like(omega)
	qdisp = np.logspace(np.log10(qmin),np.log10(qmax),num=qnum)

	qdisp,disp =PlasmonDispersionRelation(qdisp,eps1,eps2,sigma,start)

	if disp[0][0]>disp[-1][0]:
		qdisp = qdisp[::-1]
		disp = disp[::-1]

	for i in range(np.size(omega)):
		w = omega[i]

		# Eliminate omegas not in the range
		if w > np.max(disp[:,0]):
			qreturn[i] = -1
			continue
		if w < np.min(disp[:,0]):
			qreturn[i] = -1
			continue

		# Omega is in the range. Find the qdisp value closest
		qclose_i = np.argmin(np.abs(w-disp[:,0]))

		if np.abs(w-disp[qclose_i][0])/w<0.001:
			qreturn[i] = qdisp[qclose_i]
			continue

		if w > disp[qclose_i][0]:
			a = qdisp[qclose_i]
			b = qdisp[qclose_i+1]
		if omega[i] <= disp[qclose_i][0]:
			b = qdisp[qclose_i]
			a = qdisp[qclose_i-1]

		qreturn[i]=np.squeeze(InversePlasmonDispersionRelation(np.array([w]),eps1,eps2,
																sigma,start,a,b,10))

	return qreturn

def LDOSVacuum(omega):
    return omega**2 / (np.pi**2*sc.speed_of_light**3)


def RadiativeLDOSHalfspace(z,omega,eps1,eps2,sigma,field,starts=None,verbose=False):
	'''
	Returns the radiative component of the LDOS in a half space geometry.

	Currently only electric

	Parameters
	----------

	z:		height above interface (m).

	omega:	frequency at which LDOS is evaluated

	eps1:	

	sigma:	conductivity of interface. Useful with graphene.

	starts:	array of callables

	References
	----------
	[1] Joulain, Karl, Rémi Carminati, Jean-Philippe Mulet, and Jean-Jacques Greffet.
		2003. “Definition and Measurement of the Local Density of Electromagnetic States
		Close to an Interface.” Phys. Rev. B 68 (December).
		https://doi.org/10.1103/PhysRevB.68.245405.

	[2] Messina, Riccardo, Jean-Paul Hugonin, Jean-Jacques Greffet, François Marquier,
		Yannick De Wilde, Ali Belarouci, Luc Frechette, Yvon Cordier, and Philippe Ben-Abdallah. 2013.
		“Tuning the Electromagnetic Local Density of States in Graphene-Covered Systems via Strong
		Coupling with Graphene Plasmons.” Physical Review B 87 (8): 085421.
		https://doi.org/10.1103/PhysRevB.87.085421.
	'''

	if verbose:
		print('\nRadiative %s\n' %(field))

	omega = np.atleast_1d(omega)
	ldos = np.empty_like(omega)

	for i in range(np.size(omega)):
		w = omega[i]
		k0 = w/sc.speed_of_light

		rp = lambda kp: FresnelReflection(kp*k0,w,eps1,eps2,pol='p')
		rs = lambda kp: FresnelReflection(kp*k0,w,eps1,eps2,pol='s')

		# Code that finds kp at which Plasmons exist (poles)
		if starts!=None:
			for start in starts:
				radiative_poles = InversePlasmonDispersionRelation(w,eps1,eps2,sigma,
									start=start,qmin=1,qmax=k0,qnum=100)

		kperp = lambda kp: np.sqrt(1-kp**2+1e-9*1j)

		Eintegrand = lambda kp: (kp/np.abs(kperp(kp)))*( np.real(rs(kp)*np.exp(2*1j*kperp(kp)*k0*z))
	                            + np.real(rp(kp)*np.exp(2*1j*kperp(kp)*k0*z)) * (2*kp**2 - 1) )

		Hintegrand = lambda kp: (kp/np.abs(kperp(kp)))*( np.real(rp(kp)*np.exp(2*1j*kperp(kp)*k0*z))
	                            + np.real(rs(kp)*np.exp(2*1j*kperp(kp)*k0*z)) * (2*kp**2 - 1) )

		if field == 'electric':
			integrand = Eintegrand
		if field == 'magnetic':
			integrand = Hintegrand
		else:
			integrand = lambda kp: Eintegrand(kp) + Hintegrand(kp)

		integral = integrate.quad(integrand,0,1)[0]


		ldos[i] = ( LDOSVacuum(w)/4 )*(integral)

	return ldos

def EvanescentLDOSHalfspace(z,omega,eps1,eps2,sigma,field,starts=None,kcutoff=1e10,verbose=False):
	'''
	Returns the evanescent component of the LDOS in a half space geometry.

	Currently only electric.

	Parameters
	----------

	z:		height above interface (m).

	omega:	frequency at which LDOS is evaluated

	eps1:	

	sigma:	conductivity of interface. Useful with graphene.

	starts:	array of callables

	kcutoff: cutoff wavenumber for upper end of integral.

	References
	----------
	[1] Joulain, Karl, Rémi Carminati, Jean-Philippe Mulet, and Jean-Jacques Greffet.
		2003. “Definition and Measurement of the Local Density of Electromagnetic States
		Close to an Interface.” Phys. Rev. B 68 (December).
		https://doi.org/10.1103/PhysRevB.68.245405.

	[2] Messina, Riccardo, Jean-Paul Hugonin, Jean-Jacques Greffet, François Marquier,
		Yannick De Wilde, Ali Belarouci, Luc Frechette, Yvon Cordier, and Philippe Ben-Abdallah. 2013.
		“Tuning the Electromagnetic Local Density of States in Graphene-Covered Systems via Strong
		Coupling with Graphene Plasmons.” Physical Review B 87 (8): 085421.
		https://doi.org/10.1103/PhysRevB.87.085421.
	'''

	if verbose:
		print('\nEvanescent %s\n' % (field))
	omega = np.atleast_1d(omega)
	ldos = np.empty_like(omega)


	for i in range(np.size(omega)):
		w = omega[i]

		k0 = w / sc.speed_of_light
		# Define Fresnel coefficients in terms of normalized in-plane wavenumber
		rp = lambda kp: FresnelReflection(kp*k0,w,eps1,eps2,pol='p')
		rs = lambda kp: FresnelReflection(kp*k0,w,eps1,eps2,pol='s')

		# Code that finds kp at which Plasmons exist (poles)
		if starts!=None:
			for start in starts:
				poles= InversePlasmonDispersionRelation(w,eps1,eps2,sigma,
										start=start,qmin=k0,qmax=kcutoff,qnum=100)


		kperp = lambda kp: np.sqrt(1-kp**2+1e-9*1j)


		Eintegrand = lambda kp: (kp/np.abs(kperp(kp)))*( np.imag(rs(kp))
	                                                +(2*kp**2-1)*np.imag(rp(kp))) * np.exp(-2*np.abs(kperp(kp))*k0*z)

		Hintegrand = lambda kp: (kp/np.abs(kperp(kp)))*( np.imag(rp(kp))
	                                                +(2*kp**2-1)*np.imag(rs(kp))) * np.exp(-2*np.abs(kperp(kp))*k0*z)

		if field == 'electric':
			integrand = Eintegrand
		if field == 'magnetic':
			integrand = Hintegrand
		else:
			integrand = lambda kp: Eintegrand(kp) + Hintegrand(kp)

		integral = integrate.quad(integrand,1,kcutoff/k0,points=(1))[0]

		ldos[i] = ( LDOSVacuum(w)/4 )* integral

	return ldos

def LDOSHalfspace(z,omega,eps1,eps2,sigma=0,starts=None,kcutoff=1e10):
	'''

	Parameters
	----------

	z:		height above interface (m).

	omega:	frequency at which LDOS is evaluated

	eps1:	

	sigma:	conductivity of interface. Useful with graphene.

	starts:	array of callables

	kcutoff: cutoff wavenumber to end the integral.

	References
	----------
	[1] Joulain, Karl, Rémi Carminati, Jean-Philippe Mulet, and Jean-Jacques Greffet.
		2003. “Definition and Measurement of the Local Density of Electromagnetic States
		Close to an Interface.” Phys. Rev. B 68 (December).
		https://doi.org/10.1103/PhysRevB.68.245405.

	[2] Messina, Riccardo, Jean-Paul Hugonin, Jean-Jacques Greffet, François Marquier,
		Yannick De Wilde, Ali Belarouci, Luc Frechette, Yvon Cordier, and Philippe Ben-Abdallah. 2013.
		“Tuning the Electromagnetic Local Density of States in Graphene-Covered Systems via Strong
		Coupling with Graphene Plasmons.” Physical Review B 87 (8): 085421.
		https://doi.org/10.1103/PhysRevB.87.085421.
	'''

	k0 = omega/sc.speed_of_light
	rp = lambda kp: FresnelReflection(kp*k0,omega,eps1,eps2,pol='p')
	rs = lambda kp: FresnelReflection(kp*k0,omega,eps1,eps2,pol='s')

	# Code that finds kp at which Plasmons exist (poles)
	if starts:
		for start in starts:
			radiative_poles = InversePlasmonDispersionRelation(omega,eps1,eps2,sigma,
									start=start,qmin=1,qmax=k0,qnum=100)

			evanescent_poles= InversePlasmonDispersionRelation(omega,eps1,eps2,sigma,
									start=start,qmin=k0,qmax=kcutoff,qnum=100)
			print(evanescent_poles)
	kperp = lambda kp: np.sqrt(1-kp**2+1e-9*1j)

	Eintegrand1 = lambda kp: (kp/np.abs(kperp(kp)))*( np.real(rs(kp)*np.exp(2*1j*kperp(kp)*k0*z))
	                            + np.real(rp(kp)*np.exp(2*1j*kperp(kp)*k0*z)) * (2*kp**2 - 1) )
	Eintegrand2 = lambda kp: (kp/np.abs(kperp(kp)))*( np.imag(rs(kp))
	                                                +(2*kp**2-1)*np.imag(rp(kp))) * np.exp(-2*np.abs(kperp(kp))*k0*z)
	Hintegrand1 = lambda kp: (kp/np.abs(kperp(kp)))*( np.real(rp(kp)*np.exp(2*1j*kperp(kp)*k0*z))
	                            + np.real(rs(kp)*np.exp(2*1j*kperp(kp)*k0*z)) * (2*kp**2 - 1) )
	Hintegrand2 = lambda kp: (kp/np.abs(kperp(kp)))*( np.imag(rp(kp))
	                                                +(2*kp**2-1)*np.imag(rs(kp))) * np.exp(-2*np.abs(kperp(kp))*k0*z)


	Eintegral1 = integrate.quad(Eintegrand1,0,1)[0]
	Eintegral2 = integrate.quad(Eintegrand2,1,kcutoff/k0)[0]
	Hintegral1 = integrate.quad(Hintegrand1,0,1)[0]
	Hintegral2 = integrate.quad(Hintegrand2,1,kcutoff/k0)[0]

	ELDOSHalfspace = ( LDOSVacuum(omega)/4 )*(2 + Eintegral1 + Eintegral2)
	HLDOSHalfspace = ( LDOSVacuum(omega)/4 )*(2 + Hintegral1 + Hintegral2)

	LDOSHalfspace = ELDOSHalfspace + HLDOSHalfspace

	return LDOSHalfspace

def LDOS(r,omega,GE,GH):
	'''

	Eqn 9 of Ref 1


	Parameters
	---------

	r:		length 3 array, position at which to evaluate LDOS (m)

	omega:	scalar, frequency at which to evaluate LDOS (rad/s)

	GE:		DyadicGreenFunction, electric Greens function

	GH:		DyadicGreenFunction, magnetic Greens function


	References
	----------

	[1] Joulain, Karl, Rémi Carminati, Jean-Philippe Mulet, 
		and Jean-Jacques Greffet. 2003. 
		“Definition and Measurement of the Local Density of 
		Electromagnetic States Close to an Interface.” 
		Phys. Rev. B 68 (December). https://doi.org/10.1103/PhysRevB.68.245405.

	
	'''

	prefactor = omega/ (np.pi*sc.speed_of_light**2)

	tr = GE.Trace(r,r,omega) + GH.Trace(r,r,omega)

	ldos = prefactor * np.imag(tr)

	return ldos