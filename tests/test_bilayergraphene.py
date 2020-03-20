import os,sys
import pytest

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from graphenemodeling.graphene import bilayer as blg
from scipy import constants as sc
import numpy as np


class TestBilayerGraphene:

	def test_initialize(self):
		pass