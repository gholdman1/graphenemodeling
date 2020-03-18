Development
===========

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

1. All tests pass

	Change to the root directory ``graphenemodeling//`` and run ``pytest``.

2. The documentation builds

	Change to the ``doc/`` directory and run ``make``.

Linux/MacOS

.. code:: bash

	> make html

Windows

.. code:: bash

	> make.bat html

Look at the html files under the ``build`` directory.