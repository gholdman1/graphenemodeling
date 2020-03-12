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
		kF = mlg.kFermi(eF,model='LowEnergy')

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
		kF = mlg.kFermi(eF,model='LowEnergy')
		assert mlg.CarrierDispersion(kF,'LowEnergy')/sc.elementary_charge == 0.4


class TestMonolayerGraphene:

	def test_kFermi_LowEnergy(self):

		eF = 0.4 * sc.elementary_charge
		kF = mlg.kFermi(eF,model='LowEnergy')

		assert np.isclose(kF,670690811.5358821,rtol=1e-05)


	def test_DensityOfStates_LowEnergy(self):
		E = 1 * _c.g0
		DOS = mlg.DensityOfStates(E,model='LowEnergy')
		assert np.isclose(DOS,3.127898579800643e+37,rtol=1e-05)

