import os
import numpy as np

class CoupledOscillators:
	'''
	A class describing two coupled oscillators
	'''
	def __init__(self,osc1,osc2,g):
		'''
		
		Parameters
		-----------

		osc1:	length 2 array, osc1[0] is the frequency of oscillator 1,
								osc1[1] is the loss-rate of oscillator 1

		osc2:	see osc1

		g:		coupling rate
		'''

		self.osc1, self.osc2 = osc1, osc2
		self.g = g


	def Hamiltonian(self,omega):
		'''
		Hamiltonian of the coupled oscillator system.

		References
		----------

		[1] 
		'''

		H = np.array([[self.osc1[0]-omega-1j*self.osc1[1],self.g],
					  [self.g,self.osc2[0]-omega-1j*self.osc2[1]]])

		return H

	def Response(self,omega,force):
		'''
		The response of the system to an applied force f
		'''

		x = np.empty((np.size(omega),2),dtype=np.complex)

		for i in range(np.size(omega)):
			w = omega[i]
			H = self.Hamiltonian(w)
			Hinv = np.linalg.inv(H)

			x[i] = np.dot(Hinv,1j*force)


		return x