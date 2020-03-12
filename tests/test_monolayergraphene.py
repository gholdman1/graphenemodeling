import os,sys
import pytest

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import graphenemodeling.graphene.monolayer as mlg
import graphenemodeling.graphene._constants as _c
from scipy import constants as sc
import numpy as np


class TestHamiltonian:

	def test_Hamiltonian_LowEnergy(self):

		eF = 0.4 * sc.elementary_charge
		kF = mlg.FermiWavenumber(eF,model='LowEnergy')

		H = mlg.Hamiltonian(kF,model='LowEnergy')

		# np.linalg.eigh required as H is Hermitian
		energy = np.max(np.linalg.eigh(H)[0])

		assert np.isclose(energy,mlg.CarrierDispersion(kF,'LowEnergy'),rtol=1e-05)

class TestCarrierDispersion:

	def test_eh_error(self):
		with pytest.raises(ValueError):
			mlg.CarrierDispersion(0,'LowEnergy',eh=2)

	def test_CarrierDispersion_LowEnergy(self):

		eF = 0.4 * sc.elementary_charge
		kF = mlg.FermiWavenumber(eF,model='LowEnergy')
		assert mlg.CarrierDispersion(kF,'LowEnergy')/sc.elementary_charge == 0.4

class TestDensityOfStates:

	def test_DensityOfStates_LowEnergy(self):
		E = 1 * _c.g0
		DOS = mlg.DensityOfStates(E,model='LowEnergy')
		assert np.isclose(DOS,3.127898579800643e+37,rtol=1e-05)

	def test_DensityOfStates_FullTightBinding(self):
		E = 0.5 * _c.g0
		DOS = mlg.DensityOfStates(E,model='FullTightBinding')

		assert np.isclose(DOS,1.9120400733241723e+37)

	def test_DensityOfStates_FullTightBinding_inf(self):
		E = 1 * _c.g0
		DOS = mlg.DensityOfStates(E,model='FullTightBinding')

		assert DOS==np.inf

class TestkFermi:

	def test_kFermi_LowEnergy(self):

		eF = 0.4 * sc.elementary_charge
		kF = mlg.FermiWavenumber(eF,model='LowEnergy')

		assert np.isclose(kF,670690811.5358821,rtol=1e-05)

	def test_kFermi_FullTightBinding(self):

		eF = 0.4 * sc.elementary_charge
		kF = mlg.FermiWavenumber(eF,model='FullTightBinding')

		assert np.isclose(kF,671262440.2396309,rtol=1e-10)

	def test_kFermi_FullTightBinding_high_energy(self):

		eF = 0.8 * _c.g0
		kF = mlg.FermiWavenumber(eF,model='FullTightBinding')

		assert np.isclose(kF,3864007944.295662,rtol=1e-10)

class TestCarrierDensity:

	def test_LowEnergy_zero_zero(self):

		assert mlg.CarrierDensity(0,0,model='LowEnergy') == 0

	def test_LowEnergy_100_meV_0K(self):

		eV = sc.elementary_charge
		eF = 0.1 * eV

		assert np.isclose(mlg.CarrierDensity(eF,T=0,model='LowEnergy'),8949007205084713.0)

	def test_LowEnergy_zero_meV_lowT(self):

		assert np.isclose(mlg.CarrierDensity(0,T=1e-3,model='LowEnergy'),0)

	def test_LowEnergy_low_meV_lowT(self):

		eV = sc.elementary_charge
		eF = 0.01 * eV
		
		assert np.isclose(mlg.CarrierDensity(eF,T=0,model='LowEnergy'),
							mlg.CarrierDensity(eF,T=1e-3,model='LowEnergy'))

	def test_LowEnergy_100_meV_lowT(self):

		eV = sc.elementary_charge
		eF = 0.1 * eV
		
		assert np.isclose(mlg.CarrierDensity(eF,T=0,model='LowEnergy'),
							mlg.CarrierDensity(eF,T=1e-3,model='LowEnergy'))

class TestChemicalPotential:

	def test_0_0(self):
		assert mlg.ChemicalPotential(0,0)==0

	def test_1e16_0(self):
		assert np.isclose(mlg.ChemicalPotential(1e16,0),1.6936472639373315e-20)

	def test_0_lowT(self):
		assert mlg.ChemicalPotential(0,1e-3)==0

	def test_1e16_lowT(self):
		assert np.isclose(mlg.ChemicalPotential(1e16,1e-3),1.6936472639373315e-20)

	def test_1e16_10K(self):
		assert np.isclose(mlg.ChemicalPotential(1e16,10),1.693525703194968e-20)