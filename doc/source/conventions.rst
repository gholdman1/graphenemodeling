Conventions
===========

Units
-----

For ease of use and consistency, all units across the packge are SI.

This means, despite the utility of the electron-volts, units of energy of are in Joules. However, converting from eV to Joules is easy as the unit of conversion is the elementary charge.

>>> from scipy.constants import elementary_charge as eV
>>> E_eV = 0.4 # Energy in eV
>>> E_J  = E_eV * eV # Energy in Joules

Variables
---------

Chemical potentials use the variable ``mu``. If a function is only applicable at zero temperature, then we use the variable ``FermiLevel`` instead.

Examples
--------

Every piece of code must demonstrate that it replicates a published work. This is done by providing an example in the function that reproduces a plot or number from that work.

Code Structure
--------------

Just as ``numpy`` is typically imported with the name ``np``, we have a few conventions in GrapheneModeling. The submodules ``graphene.monolayer`` and ``graphene.bilayer`` are imported as ``mlg`` and ``blg``.

>>> from graphenemodeling.graphene import monolayer as mlg
>>> from graphenemodeling.graphene import bilayer as blg

One logical use case of this code is to analyze a piece of graphene at a temperature :math:`T` with mobility :math:`\\mu`. In this case, the submodule ``graphene.monolayer`` sould have been a class, perhaps called ``Monolayer`` with attributes ``Monolayer.T`` and ``Monolayer.mobility``. However, since most experiments involve fitting to the properties of a piece of graphene with unknown mobility and temperature, this approach would have been slow, as the user would be required to reinstantiate the ``Monolayer`` class for every mobility and temperature of interest.
