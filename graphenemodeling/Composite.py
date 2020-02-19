import os
import numpy as np
import warnings

import fundamental_constants as fc
import statistical_distributions as sd
import Electrodynamics as em

class Stack:

	def __init__(self,layers,T):
		'''

		Parameters
		----------

		background:	length 2 array. Semi-infinite background on
					either side of dielectric.

		layers:		A list of finite-sized layers
					Format: [ [material1, thickness1],
							  [material2, thickness2],
							  ...,
							  [materialn,thicknessn]]

		T:			Temperature (K)
		
		'''

		self.Temperature = T
		self.total_thickness = 0

		if layers[0][1] != np.inf:
			raise Exception('Must specify material for z<0. Usage: layers[0]=[material, np.inf]')

		if layers[-1][1] != np.inf:
			raise Exception('Must specify material above. Usage: layers[-1]=[material, np.inf]')


		self.layer=[]

		for i in range(len(layers)):
			if i !=0 and i !=len(layers)-1:
				if layers[i][1] == np.inf:
					raise Exception('Infinite layer at Layer %d' % (i) )

				self.total_thickness=self.total_thickness + layers[i][1]

			self.add_layer(layers[i][0],
						   layers[i][1])

	def ComputeReflectance(self):
		'''
		The reflectance of the dielectric stack.
		'''

		r, t = em.FresnelCoefficientStack()