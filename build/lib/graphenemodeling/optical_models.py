import os
import numpy as np
from scipy import integrate

def LorentzTerm(w,w0,gamma,s):
	'''
	A single term in a Lorentz Oscillator sum.

	Parameters
	----------

	w:		array-like, frequency at which to evaluate the response (rad/s)

	w0:		scalar, resonance frequency (rad/s)

	gamma:	scalar, width of the resonance (rad/s)

	s:		scalar, strength of the response (unitless)
	'''

	res = s * w0**2 / (w0**2 - w**2 - 1j*gamma*w)

	return res

def Lorentzian(omega,epsinf,modes):
	'''
	The full Lorentzian model.

	Parameters
	----------

	epsinf:		scalar, the infinite frequency permittivity

	modes:		array-like, n x 3 array of responses
							Format:[[w01,gamma1,s1],
									[w02,gamma2,s2],
									...]
	'''

	LorentzSum = np.zeros_like(omega)

	for i in range(np.shape(modes)[0]):

		LorentzSum = LorentzSum + LorentzTerm(omega,modes[i][0],
													modes[i][1],
													modes[i][2])

	return epsinf + LorentzSum

def DrudeSommerfeld(omega,omegap,gamma):
	'''
	The Drude-Sommerfeld model of a metal dielectric.
	Applies when only free-electron contributions need
	to be considered.

	Returns
	----------

	eps:		Drude permittivity
	'''

	frac = omegap**2 / (omega**2 + 1j*gamma*omega)

	eps = 1 - frac

	return eps

def Fano(omega,omega0,gamma,q):
	'''
	A Fano lineshape

	Parameters
	------------

	omega:		array-like, frequency (rad/s)

	omega0:		scalar, resonant frequency (rad/s)

	gamma:		scalar, loss rate (rad/s)
	
	q:			scalar, Fano asymmetry parameters

	References
	---------

	[1] http://demonstrations.wolfram.com/FanoResonance/.

	'''

	wReduced = 2 * (omega-omega0)/gamma

	numerator = (wReduced+q)**2

	denominator = wReduced**2+1

	fano = numerator / denominator

	return fano

def Brendel(omega,omega0,sd0,omegap,omegat):
	'''
	A broadened version of the Lorentz model.

	Equation 3 of Ref 1

	Parameters
	----------

	Returns
	----------


	References
	----------

	[1] Kischkat, Jan, Sven Peters, Bernd Gruska, Mykhaylo Semtsiv,
		Mikaela Chashnikova, Matthias Klinkmüller, Oliana Fedosenko,
		et al. 2012.
		“Mid-Infrared Optical Properties of Thin Films of Aluminum Oxide,
		Titanium Dioxide, Silicon Dioxide, Aluminum Nitride,
		and Silicon Nitride.” Applied Optics 51 (28): 6789–98.
		https://doi.org/10.1364/AO.51.006789.

	'''

	Gaussian = lambda x, x0, sd: np.exp((x-x0)**2/(2*sd**2)) / (sd*np.sqrt(2*fc.pi))

	integrand= lambda x,omega: Gaussian(x,omega0,sd0) * omegap**2 / (x**2 - omega**2-1j*omegat*omega)

	integral = integrate.quad(integrand,-np.inf,np.inf,args=(omega))

	return integral
