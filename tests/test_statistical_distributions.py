import os,sys
import pytest

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import graphenemodeling.statistical_distributions as sd
from scipy import constants as sc


class TestFermiDirac:

	def test_zero(self):
		fd = sd.FermiDirac(0,10)

		assert fd == 0.5

		