Introduction
============

Purpose
-------

The GrapheneModeling package is a collection of Python code intended to model measurements of graphene.

Getting Started
---------------

To use GrapheneModeling, you may pip install it.

.. code:: bash

	> pip install graphenemodeling

Since GrapheneModeling is dependent on Numpy, SciPy, and Matplotlib, this will install these packages as well.
You may also want to install Jupyter if you prefer working in that sort of environment.

Test out this installation by running the following

.. code:: bash

	> python -m graphenemodeling

You should see a welcome message

Submodules
----------

GrapheneModeling is organized into a few submodules. The first is ``graphene``. 
This contains the theoretical expressions for the properties of graphene.

Just as ``numpy`` is typically imported with the name ``np``, we have a few conventions in GrapheneModeling. The submodules ``graphene.monolayer`` and ``graphene.bilayer`` are imported as ``mlg`` and ``blg``.

>>> from graphenemodeling.graphene import monolayer as mlg
>>> from graphenemodeling.graphene import bilayer as blg

Units
-----

For ease of use and consistency, all units across the package are SI. The units of parameters are also documented in every function.

This means, despite the utility of the electron-volt (eV), units of energy of are in Joules (J). Converting from eV to J is easy, however, as the unit of conversion is the elementary charge :math:`e`.

>>> from scipy.constants import elementary_charge as eV
>>> E_eV = 0.4 # Energy in eV
>>> E_J  = E_eV * eV # Energy in Joules

Variable Names
--------------

Chemical potentials use the variable ``mu``. If a function is only applicable at zero temperature, then we use the variable ``FermiLevel`` instead. This latter variable is only used when there is no ambiguity, like in the function ``FermiWavenumber``, which only exists when temperature is zero.

The variable ``k`` is used for the wavenumber of carriers such as Dirac fermions. The variable ``q`` is used for scattering wavevectors and the quantities derived from such scattering considerations such as ``Polarizibility``, ``OpticalConductivity`` and ``PlasmonDispersion``.

Examples
--------

Every function in GrapheneModeling should be accompanied by an example replicating a result in the literature.
This allows users to a) be assured the code is accurately written, b) understand how to use the function, and c)
have a reference on where to find more information.