import numpy as np

from scipy import integrate


class Ellipsoid:

	def __init__(self,eps,semimajor):

		self.eps = eps
		self.ax = semimajor[0]
		self.ay = semimajor[1]
		self.az = semimajor[2]

	def ShapeFactor(self):
		'''
		
		Equation 6 of Ref 1

		References
		-----------

		[1] Sihvola, A. (2007). Dielectric Polarization and Particle Shape Effects.

		'''
		ax, ay, az = self.ax, self.ay, self.az

		prefactor = ax * ay * az / 2

		sfactor = lambda s, a: np.sqrt( s + a**2 )

		integrand = lambda s, ax,ay,az: 1 / ( sfactor(s,ax)**3 * sfactor(s,ay) * sfactor(s,az) )

		Nx = prefactor*integrate.quad(integrand,0,np.inf,args=(ax,ay,az))[0]

		return Nx

	def Polarizibility(self,MLWA=False,):
		'''
		
		Equation 8 of Ref 1

		Parameters
		----------

		MLWA:		Bool, whether to use the modified long-wavelength approximation.

		References
		-----------

		[1] Sihvola, A. (2007). Dielectric Polarization and Particle Shape Effects.

		'''

		Nx = self.ShapeFactor()

		alpha = (self.eps - 1) / ( 1 + Nx*(self.eps-1) )
		

		if MLWA:
			'''
			Code not written
			'''
			pass

		return alpha


