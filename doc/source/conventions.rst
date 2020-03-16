Conventions
===========

Units
-----

For ease of use and consistency, all units across the package are SI. The units of parameters are also documented in every function.

This means, despite the utility of the electron-volt (eV), units of energy of are in Joules (J). Converting from eV to J is easy, however, as the unit of conversion is the elementary charge :math:`e`.

>>> from scipy.constants import elementary_charge as eV
>>> E_eV = 0.4 # Energy in eV
>>> E_J  = E_eV * eV # Energy in Joules

Variables
---------

Chemical potentials use the variable ``mu``. If a function is only applicable at zero temperature, then we use the variable ``FermiLevel`` instead. This latter variable is only used when there is no ambiguity, like in the function ``FermiWavenumber``, which only exists when temperature is zero.

The variable ``k`` is used for the wavenumber of carriers such as Dirac fermions. The variable ``q`` is used for scattering wavevectors and the quantities derived from such scattering considerations such as ``Polarizibility``, ``OpticalConductivity`` and ``PlasmonDispersion``.

Examples
--------

Every piece of code must demonstrate that it replicates a published work. This is done by providing an example in the function that reproduces a figure or number from that work.

Note also, that ``numpy`` is implicitly imported in examples. So when you replicate the example, you will need to add a line

.. code:: python

	>>> import numpy as np


Code Structure
--------------

Just as ``numpy`` is typically imported with the name ``np``, we have a few conventions in GrapheneModeling. The submodules ``graphene.monolayer`` and ``graphene.bilayer`` are imported as ``mlg`` and ``blg``.

>>> from graphenemodeling.graphene import monolayer as mlg
>>> from graphenemodeling.graphene import bilayer as blg

One logical use case of this code is to analyze a piece of graphene at a temperature :math:`T` with mobility :math:`\mu`. In this case, the submodule ``graphene.monolayer`` sould have been a class, perhaps called ``Monolayer`` with attributes ``Monolayer.T`` and ``Monolayer.mobility``. However, since most experiments involve fitting to the properties of a piece of graphene with unknown mobility and temperature, this approach would have been slow, as the user would be required to reinstantiate the ``Monolayer`` class for every :math:`\mu` and :math:`T` of interest.
