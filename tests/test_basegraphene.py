import os,sys
import pytest

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from graphenemodeling.graphene import BaseGraphene
from scipy import constants as sc


class TestBaseGraphene:

	def test_initialize(self):
		bg = BaseGraphene()

	def test_constants(self):
		bg = BaseGraphene()

		assert bg.a == 1.42e-10
		assert bg.A == 3*(3**(1/2))*(bg.a**2) / 2
		bg.g0 == 2.8*sc.elementary_charge
		bg.g0prime == 0 * bg.g0
		bg.vF == 3*bg.a*bg.g0/(2*sc.hbar)