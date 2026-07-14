Installation
============

CompactObject is an open-source package designed for comprehensive neutron star EOS inference. It is built to be easy to install and use. Follow the step-by-step installation guide below to get started.

Using Python Virtual Environment (Recommended)
----------------------------------------------

If you are not using Anaconda, you can create a virtual environment using Python's ``venv`` module:

Step 1: Create a Virtual Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the following command to create a virtual environment named ``CompactObject``:

.. code-block:: bash

   python3 -m venv CompactObject

*You can specify a different path by replacing ``CompactObject`` with your desired directory name.*

Step 2: Activate the Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Activate the virtual environment with:

.. code-block:: bash

   source CompactObject/bin/activate

Step 3: Install the CompactObject Package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the environment is activated, you can install the package from PyPI, which is called ``CompactObject-TOV``:

.. code-block:: bash

   pip install CompactObject-TOV

The dependencies should be automatically installed for you.
This includes the packages used by the shipped examples, such as ``pandas``,
``sympy``, ``ultranest``, ``h5py``, and ``numba``. The only accelerated
dependency that is intentionally optional is ``NumbaMinpack``; see
`Optional FastRMF dependency`_ below.

Alternative Installation Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you encounter issues using ``pip install``, you can install the package manually:

1. **Clone the Repository**

   Clone the CompactObject repository:

   .. code-block:: bash

      git clone https://github.com/ChunHuangPhy/CompactObject.git

2. **Install the Required Dependencies**

   Navigate to the repository directory and install the dependencies:

   .. code-block:: bash

      pip install -r requirements.txt

3. **Install the CompactObject Package**

   Install the package in editable mode:

   .. code-block:: bash

      pip install -e .

   To run the example notebooks from a source checkout, install the notebook
   extras as well. The base package supplies the scientific runtime
   dependencies; this extra supplies Jupyter and notebook execution tools:

   .. code-block:: bash

      pip install -e ".[notebooks]"

   To include the optional FastRMF benchmark path used by some examples, use:

   .. code-block:: bash

      pip install -e ".[notebooks,fast]"

To upgrade to the latest version on PyPI, run:

.. code-block:: bash

   pip install CompactObject-TOV --upgrade

Step 4: Using the Package
^^^^^^^^^^^^^^^^^^^^^^^^^

You are now ready to use CompactObject. Each time you want to use the package, ensure you activate the environment:

.. code-block:: bash

   source CompactObject/bin/activate

Alternative: Using Anaconda
---------------------------

Step 1: Create a Virtual Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also create a virtual environment for CompactObject using Anaconda:

.. code-block:: bash

   conda create -n CompactObject

When prompted to proceed, type ``y`` and press Enter.

Step 2: Activate the Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Activate your newly created environment with the following command:

.. code-block:: bash

   conda activate CompactObject

**Note:** Once you create this environment, you don't need to create it again. Simply activate it whenever you want to use CompactObject.

Step 3: Install the CompactObject Package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the environment is activated, you can install the package from PyPI:

.. code-block:: bash

   pip install CompactObject-TOV

The dependencies should be automatically installed for you.
This includes the packages used by the shipped examples, such as ``pandas``,
``sympy``, ``ultranest``, ``h5py``, and ``numba``. The only accelerated
dependency that is intentionally optional is ``NumbaMinpack``; see
`Optional FastRMF dependency`_ below.

Alternative Installation Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you encounter issues using ``pip install``, you can install the package manually:

1. **Clone the Repository**

   Clone the CompactObject repository:

   .. code-block:: bash

      git clone https://github.com/ChunHuangPhy/CompactObject.git

2. **Install the Required Dependencies**

   Navigate to the repository directory and install the dependencies:

   .. code-block:: bash

      pip install -r requirements.txt

3. **Install the CompactObject Package**

   Install the package in editable mode:

   .. code-block:: bash

      pip install -e .

   To run the example notebooks from a source checkout, install the notebook
   extras as well. The base package supplies the scientific runtime
   dependencies; this extra supplies Jupyter and notebook execution tools:

   .. code-block:: bash

      pip install -e ".[notebooks]"

   To include the optional FastRMF benchmark path used by some examples, use:

   .. code-block:: bash

      pip install -e ".[notebooks,fast]"

To upgrade to the latest version on PyPI, run:

.. code-block:: bash

   pip install CompactObject-TOV --upgrade

Step 4: You're Ready to Use CompactObject!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Whenever you want to use this package, remember to activate the environment first:

.. code-block:: bash

   conda activate CompactObject

By following these instructions, you should have CompactObject installed and ready to use. If you encounter any issues during installation, please refer to the project's documentation or seek assistance from the community.

..    Our package automatically installs all necessary dependencies for you. The dependencies include:

..    - `corner`
..    - `csv`
..    - `itertools`
..    - `math`
..    - `matplotlib`
..    - `numba`
..    - `numbaminpack`
..    - `numpy`
..    - `os`
..    - `pandas`
..    - `scipy`
..    - `sys`
..    - `ultranest`

Optional FastRMF dependency
---------------------------

``CompactObject-TOV`` can use ``NumbaMinpack`` for the accelerated
``EOSgenerators.fastRMF_EoS`` path. This dependency is optional; the standard
RMF and DDH examples work without it.

``NumbaMinpack`` requires a Fortran compiler before installation. On macOS,
the upstream project documents the following setup:

.. code-block:: bash

   brew install gcc
   python -m pip install NumbaMinpack

On Debian/Ubuntu systems, install ``gfortran`` first:

.. code-block:: bash

   sudo apt-get install gfortran cmake
   python -m pip install NumbaMinpack

To install CompactObject with the optional FastRMF dependencies from a source
checkout, run:

.. code-block:: bash

   pip install -e ".[fast]"

To install everything needed for documentation and notebook execution:

.. code-block:: bash

   pip install -e ".[docs,notebooks,fast]"

.. Summary
.. -------

.. - **Using Anaconda:**
..   1. Create and activate the `CompactObject` environment.
..   2. Install CompactObject with `pip`.
..   3. Activate the environment whenever you use the package.

.. - **Using Python Virtual Environment:**
..   1. Create and activate the `CompactObject` virtual environment.
..   2. Install CompactObject with `pip`.
..   3. Activate the environment whenever you use the package.

If you encounter any issues or have questions, feel free to reach out for support. Happy computing!
