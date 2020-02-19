import numpy as np

# Script that defines fundamental constants
# All constants are in SI units

e_proton		= 1.602176*1e-19		# C 		, Fundamental charge
hbar			= 1.05457148*1e-34		# J-s		, Reduced Planck's constant
kB				= 1.3806503*1e-23		# J/K 		, Boltzmann Constant
epsilon_0		= 8.854187817*1e-12 	# F/M		, Vacuum permittivity
eps0 			= 8.854187817*1e-12		# F/M 		, Vacuum permittivity
mu0				= 1.2566370614*1e-6 	# Wb/(A-m)	, Vacuum permeability
alpha			= (137.035999)**(-1) 	# 			, Fine structure constant
pi 				= 3.1415926535			# 			, pi

# Derived quantities
Z0				= (mu0/epsilon_0)**(1/2)	#  		, Vacuum impedance
c0				= (mu0*epsilon_0)**(-1/2)  # 		, speed of light in vacuum

sigma_0			= e_proton**2 / (4*hbar)	#		, Fundamental conductivity (2D?)