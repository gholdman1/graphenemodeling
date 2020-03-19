Development
===========

If you have expertise or code that fits the Purpose Statement, then you are encourage to help develop GrapheneModeling.

Development takes place on `Github <https://github.com/gholdman1/graphenemodeling/>`_.

If you are new to development, you may find the page :doc:`development.from_zero` to be useful.


Getting Started
---------------

Clone the repository

``> git clone https://gitub.com/gholdman1/graphenemodeling.git .``

Change directories into the repository.

``> cd graphenemodeling``

Development takes place on the ``develop`` branch.

``> git checkout develop``

Make a branch off of this and submit changes.

Commit Messages
---------------

This project uses the same commit acronyms as SciPy. The relevant ones are listed here.

.. code:: bash

	API: an (incompatible) API Change
	BUG: bug fix
	DEP: deprecate something, or remove a deprecated object
	DOC: documentation
	ENH: enhancement
	MAINT: maintenance commit (refactoring, typos, etc.)
	REV: revert and earlier commit
	STY: style fix (whitespace, PEP8)
	TST: addition or modification of tests
	REL: related to releasing GrapheneModeling

Checklist
---------

Before submitting ensure

1. If you added a new function:
	
	a. Your parameters are documented with units. They should be SI.

	b. The return value and units are documented.

	c. An example is included replicating a piece of literature.

	d. References are included.

	e. Possible errors in usage are caught by raising exceptions

2. Tests are included and pass
	
	Tests are housed in the ``tests/`` directory. Add class for testing your function
	and test that values are returned correctly and, if applicable, exceptions are caught.

	Change to the root directory ``graphenemodeling//`` and run ``pytest``.

2. The documentation builds

	Change to the ``doc/`` directory and run ``make html`` (Linux/MacOS) or ``make.bat html`` (Windows)

	Look at the html files under the ``build`` directory. Ensure especially that your example is plotted.