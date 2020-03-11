import numpy as np
import warnings

import scipy.constants as sc 

class Dipole:
	'''

	References
	----------

	[1] Novotny, L., and Hecht, B. (2006). Principles of Nano-Optics (Cambridge University Press).
		http://www.fulviofrisone.com/attachments/article/406/Principles%20of%20Nano-Optics.pdf

	'''
	def __init__(self,d,q0):

		self.Dipole=d
		self.IntrinsicQuantumYield=q0
		pass

	def FreeSpaceDecayRate(self,omega,d):
		'''
		Free space decay rate of a dipole emitter

		Equation 8.120 of Ref 1.

		'''
		warnings.warn('Dipole.DecayRate has not been validated')

		return (omega/sc.speed_of_light)**3 * np.linalg.norm(d,d)**2 / (3*np.pi*sc.epsilon_0*sc.hbar)

	def FreeSpacePower(self,omega,d):
		'''
		Power radiated by dipole emitter

		Equation 8.126 of Ref 1
		'''

		return self.FreeSpaceDecayRate*sc.hbar*omega
