import os
import numpy as np
import warnings

import fundamental_constants as fc
import optical_models as om


class Metal:

	'''
	Outline for a class describing metals.

	'''

	def __init__(self, model,modelparams):

		if model=='DrudeSommerfeld':
			if np.size(modelparams)!=2:
				raise Exception('modelparams must have exactly two values:\n np.array([plasma frequency, loss rate]')
			self.PlasmaFrequency = modelparams[0]
			self.Gamma  = modelparams[1]
			self.Permittivity = lambda kp, w: om.DrudeSommerfeld(w,self.PlasmaFrequency,
																	self.Gamma)

		else:
			raise Exception('Permittivity model %s is unavailable' % (model))
class Dielectric:

	def __init__(self):
		pass


	def Permittivity(self,kpar,omega):
		pass


class Gold(Metal):

	def __init__(self):

		self.modelparams = np.array([13.8e15, 1.075e14])
		Metal.__init__(self,'DrudeSommerfeld',self.modelparams)

class Aluminum(Metal):

	def __init__(self):
		'''
		
		References
		----------

		[1] Palik
		'''
		self.modelparams = np.array([1.747e16,7.596e13])

		Metal.__init__(self,'DrudeSommerfeld',self.modelparams)

class SiliconCarbide:

	def __init__(self):
		'''
		From Spitzer et al. oscillator model

		omegaL: 969 cm-1 = 1.827e14 rad/s

		omegaT: 793 cm-1 = 1.495e14 rad/s

		Gamma: 4.76 cm-1 = 0.9e12 rad/s

		'''
		self.epsinf = 6.7
		self.modelparams = np.array([1.827e14,1.495e14,0.9e12])

		self.wspp = 1.787e14 # Surface plasma frequency


	def Permittivity(self,q,omega):
		'''
		Permittivity of SiC as given by Spitzer et al.
		'''
		num = ( self.modelparams[0]**2 - self.modelparams[1]**2 )

		den = self.modelparams[1]**2-omega**2-1j*self.modelparams[2]*omega

		eps = self.epsinf * (1 + num/den)

		return eps
		
class HexagonalBoronNitride(Dielectric):

	def __init__(self,model):
		'''
		

		References
		----------

		[1] Geick, R., C. H. Perry, and G. Rupprecht. 1966.
			“Normal Modes in Hexagonal Boron Nitride.”
			Physical Review 146 (2): 543–47.
			https://doi.org/10.1103/PhysRev.146.543.

		[2] Woessner, Achim, Mark B. Lundeberg, Yuanda Gao,
			Alessandro Principi, Pablo Alonso-González,
			Matteo Carrega, Kenji Watanabe, et al. 2015.
			“Highly Confined Low-Loss Plasmons in Graphene–Boron
			Nitride Heterostructures.” Nature Materials 14 (4): 421–25.
			https://doi.org/10.1038/nmat4169.

		[3] Cai, Yongqing, Litong Zhang, Qingfeng Zeng, Laifei Cheng,
			and Yongdong Xu. 2007. “Infrared Reflectance Spectrum of BN
			Calculated from First Principles.”
			Solid State Communications 141 (5): 262–66.
			https://doi.org/10.1016/j.ssc.2006.10.040.

		[3] Brar, Victor W., Min Seok Jang, Michelle Sherrott, Seyoon Kim,
			Josue J. Lopez, Laura B. Kim, Mansoo Choi, and Harry Atwater. 2014.
			“Hybrid Surface-Phonon-Plasmon Polariton Modes in Graphene/Monolayer
			h-BN Heterostructures.” Nano Letters 14 (7): 3876–80.
			https://doi.org/10.1021/nl501096s.

		'''

		if model=='Cai':
			pass

		if model=='Cai:clean':
			self.epsinf_xy 	= 4.87
			self.epsinf_z	= 2.95

			# Strength
			self.s_xy	= [1.83]
			self.s_z	= [0.61]

			# Frequency
			self.w_xy	= [0.1701*fc.e_proton /fc.hbar]
			self.w_z	= [0.0925*fc.e_proton /fc.hbar]

			# Loss
			self.g_xy	= [0.00087*fc.e_proton/fc.hbar]
			self.g_z	= [0.00025*fc.e_proton/fc.hbar]

	def PermittivityInPlane(self,omega):

		modes = [[self.w_xy[0],self.g_xy[0],self.s_xy[0]]]
		eps = om.Lorentzian(omega,self.epsinf_xy,modes)

		return eps

	def PermittivityOutOfPlane(self,omega):
		pass

	def Permittivity(self,omega):
		epsx = self.PermittivityInPlane(omega)
		epsy = epsx
		epsz = self.PermittivityOutOfPlane(omega)

		return np.diag(epsx,epsy,epsz)

####################
# Useful Functions #
####################

def download_material_data(url,material,filename):
	'''
	Download data from a website, i.e. refractiveindex.info
	'''

	from urllib import request

	savepath = os.path.join(os.environ['DATA'],'materials',material,filename)
	request.urlretrieve(url,savepath)

def get_material_data_files(material):
	'''

	'''

	path = os.path.join(os.environ['DATA'],'materials',material)

	return os.listdir(path)

def load_material_data(material,filename):
	'''
	Loads a CSV file of data to a numpy array
	'''

	path = os.path.join(os.environ['DATA'],'materials',material,filename)

	with open(path) as f:

		data = np.loadtxt(f,delimiter=',',skiprows=1)

	return data
