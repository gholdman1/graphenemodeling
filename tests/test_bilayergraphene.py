import os,sys
import pytest

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from graphenemodeling.graphene.bilayer import Bilayer
from scipy import constants as sc
import numpy as np


class TestBilayerGraphene:

	def test_initialize(self):
		blg = Bilayer()