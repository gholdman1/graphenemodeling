import os,sys
import pytest

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import graphenemodeling.statistical_distributions as sd
from scipy import constants as sc
import numpy as np


class TestFermiDirac:

	def test_zero(self):
		fd = sd.FermiDirac(0,10)

		assert fd == 0.5

		
class TestBoseEinstein:

	def test_single_value(self):
		kB = sc.k
		T=300
		E = kB*T


		#be = sd.BoseEinstein(E,T)

