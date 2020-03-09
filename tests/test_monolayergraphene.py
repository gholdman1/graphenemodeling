import os,sys
import pytest

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from graphenemodeling.graphene import Monolayer
from scipy import constants as sc
import numpy as np


class TestMonolayerGraphene:

	def test_kFermi_LowEnergy(self):
		mlg = Monolayer()

		eF = 0.4 * sc.elementary_charge
		kF = mlg.kFermi(eF,model='LowEnergy')

		assert np.isclose(kF,670690811.5358821)

	def test_DiracFermionDispersion_LowEnergy(self):
		mlg = Monolayer()

		eF = 0.4 * sc.elementary_charge
		kF = mlg.kFermi(eF,model='LowEnergy')
		assert mlg.DiracFermionDispersion(kF,'LowEnergy')/sc.elementary_charge == 0.4

	def test_DensityOfStates_LowEnergy(self):
		mlg = Monolayer()
		E = 1 * mlg.g0
		DOS = mlg.DensityOfStates(E,model='LowEnergy')
		assert np.isclose(DOS,3.127898579800643e+37)