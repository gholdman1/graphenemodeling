Development from Zero
=====================

This is page is for individuals who have never contributed to an open-source project before.

Make a Github Account
---------------------

Go to `Github <https://github.com/join?source=header-home/>`_ and make an account.

Downloading Development Software
--------------------------------

You will need some specific pieces of software before starting.

**1) Git**

Git is a very common command line tool used to keep track of changes to projects.
It is how you will receive this project and tell the maintainer about the changes you made.

On **MacOS** and **Linux**, it should be easily available/installable from the Terminal.

On **Windows**, you will need to download `Git Bash <https://git-scm.com/downloads/>`_. This is a terminal that contains most of the 
familiar bash commands (``cd``, ``rm``, etc.) along with git itself.

**2) Conda**

You will also need to download a package manager. I recommend ``conda``. This is easily installable from a terminal on **Linux**. If you are on **Windows**, Download miniconda_ here. This will also install a terminal with the ``conda`` package manager installed.

.. _miniconda: https://docs.conda.io/en/latest/miniconda.html

Obtaining the Project
---------------------

Sign in to Github and go to `this project's page <https://github.com/gholdman1/graphenemodeling/>`_.

Click the ``Fork`` button at the top right (next to ``Star``). This will create a copy on your account.

Open your terminal containing ``git``.

Make a directory that's a sensible location for this project. For example, ``Programming\OpenSource\``.

Change into this directory.

``> cd Programming\OpenSource\``

Make sure you are connected to the internet. Run the following command.

``> git clone https://github.com/<your-github-username>/graphenemodeling.git``

This will download your version of the project. You should now see a directory ``graphenemodeling``.

Creating a Development Environment
----------------------------------

Open ``conda``. If you downloaded the miniconda installer, open the ``Anaconda Prompt`` that it installed on your computer.

Create a new environment and install the packages required for ``graphenemodeling`` to run.

.. code:: bash

	> conda create -n graphenemodeling-dev python=3.7 scipy matplotlib

Now activate this environment

.. code:: bash

	> conda activate graphenemodeling-dev

Install ``pytest`` for testing and ``sphinx`` for building the documentation.

.. code:: bash

	> conda install pytest sphinx

Finally, we need to tell the environment where ``graphenemodeling`` lives. Change directories to the root directory of the project.

.. code:: bash

	> cd Programming\OpenSource\graphenemodeling

And run

.. code:: bash

	> python setup.py develop

Run the following to ensure everything works.

.. code:: bash

	> python -m graphenemodeling


Making changes to the Project
-----------------------------

When developing, I like to have two terminals open: a testing terminal (``conda``) and a versioning terminal (``git``). Open these two.

Git has a feature called "branches", which are essentially just different versions of the code.
The released version of the project (the one on PyPI) is on the ``master`` branch. You should not make changes to this branch as they will not be accepted.

The developing version of the project is on the ``develop``
branch. All changes you make should start from the ``develop`` branch.
To accomplish this, checkout the ``develop`` branch.

.. code:: bash

	> git checkout develop


Then, make a new branch and check it out.

.. code:: bash

	> git branch new-feature
	> git checkout new-feature


As long as the ``new-feature`` branch is checked out, you can be bold with your changes.

Any time you start a new session working on the project, follow this workflow.

.. code:: bash

	# Ensure the master branch and development branch are up-to-date
	> git checkout master
	> git fetch upstream
	> git merge
	> git checkout develop
	> git merge
	# Checkout your feature branch
	> git checkout new-feature
	# Rebase so that all changes replay as if they're on the head of develop
	> git rebase develop

Then repeat the following workflow as many times as you need.

.. code:: bash

	# Make more changes (typeity type type)
	# Check that the files you wanted to change have changed
	> git status
	# Add the files
	> git add .
	> git commit -m "Type a message here. What did you change?"

Finally, when you are ready to submit your changes,

.. code:: bash

	# Push the new-feature branch to your account
	> git push

Go to your github account and make a pull request.