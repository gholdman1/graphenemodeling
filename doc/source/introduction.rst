Introduction
============

Purpose
-------

The GrapheneModeling package is a collection of Python code modeling propeties and measurements of graphene built on references to graphene literature. By putting cited descriptions of graphene's properties into the docstrings of Python functions, we hope to

1. Create a graphene wiki: a quick reference for inviduals discussing graphene's properties.
2. Clarify models used to calculate graphene properties and measurements.
3. Reduce ramp-up time for researchers requiring code to model graphene.

Getting Started
---------------

To use GrapheneModeling, you may pip install it.

.. code:: bash

	> pip install graphenemodeling

Since GrapheneModeling is dependent on Numpy, SciPy, and Matplotlib, this will install or update these packages as well. Use an environment if you'd like to keep the packages contained.
You may also want to install Jupyter if you prefer working in a notebook environment.

Test out this installation by running the following from a command line.

.. code:: bash

	> python -m graphenemodeling

You should see a welcome message.

Units
-----

For ease of use and consistency, all units across the package are SI. The units of parameters are also documented in every function.

This means, despite the utility of the electron-volt (eV), units of energy of are in Joules (J). Converting from eV to J is easy, however, as the unit of conversion is the elementary charge :math:`e`.

>>> from scipy.constants import elementary_charge as eV
>>> E_eV = 0.4 # Energy in eV
>>> E_J  = E_eV * eV # Energy in Joules

Submodules
----------

GrapheneModeling is organized into submodules. You may view them at the :doc:`documentation` page.

Just as ``numpy`` is typically imported with the name ``np``, we have a few conventions in GrapheneModeling. The submodules ``graphene.monolayer`` and ``graphene.bilayer`` are imported as ``mlg`` and ``blg``.

>>> from graphenemodeling.graphene import monolayer as mlg
>>> from graphenemodeling.graphene import bilayer as blg

Variable Names
--------------

For clarity, GrapheneModeling reserves certain variable names for certain quantities. Please raise an `issue on Github <https://github.com/gholdman1/graphenemodeling/issues/>`_ if you find any code is inconsistent with the rules outlined here.

Chemical potentials use the variable ``mu``. If a function is only applicable at zero temperature, then we use the variable ``FermiLevel`` instead. This latter variable is only used when there is no ambiguity, like in the function ``FermiWavenumber``, which only exists when temperature is zero.

The variable ``k`` is used for the wavenumber of carriers such as Dirac fermions. The variable ``q`` is used for scattering wavevectors and the quantities derived from such scattering considerations such as ``Polarizibility``, ``OpticalConductivity`` and ``PlasmonDispersion``.

Permittivities :math:`\\epsilon` are denoted by ``eps`` and will all have units of F/m. Dielectric constants (aka "relative permittivities") :math:`\\kappa` are denoted by ``kappa`` and are unitless.

Examples
--------

Every function in GrapheneModeling should be accompanied by an example replicating a result in the literature. This allows users to

1. Be assured the code is accurately written
2. Understand how to use the function
3. Have a reference on where to find more information.

In order to learn how to use GrapheneModeling, follow along with the examples in the documentation.