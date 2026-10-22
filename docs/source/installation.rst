Installation
============

CompactObject supports installation with either Python's ``venv`` module or
Conda. The package name on PyPI is ``CompactObject-TOV``.

Python virtual environment
--------------------------

Create and activate an isolated environment:

.. code-block:: bash

   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip

On Windows, activate the environment with ``.venv\Scripts\activate`` instead.

Install the released package:

.. code-block:: bash

   python -m pip install CompactObject-TOV

Release `2.1.1 <https://pypi.org/project/CompactObject-TOV/2.1.1/>`_ includes
the standard example-notebook dependencies, including ``pandas``, ``sympy``,
``ultranest``, ``h5py``, and ``numba``. If CompactObject is already installed,
upgrade it in the same Python environment before running the notebooks:

.. code-block:: bash

   python -m pip install --upgrade CompactObject-TOV

The optional FastRMF dependency ``NumbaMinpack`` is installed separately; see
:ref:`optional-fastrmf-installation` for its prerequisites and GitHub
installation instructions.

To run the example notebooks from a source checkout, install the package's
``notebooks`` extra:

.. code-block:: bash

   git clone https://github.com/ChunHuangPhy/CompactObject.git
   cd CompactObject
   python -m pip install -e ".[notebooks]"

The dependencies declared in ``pyproject.toml`` are installed automatically.
There is no separate ``requirements.txt`` dependency list.

Conda environment
-----------------

The repository includes an ``environment.yml`` file. It explicitly installs
Python into the environment and installs the runtime and notebook dependencies
from ``conda-forge``:

.. code-block:: bash

   git clone https://github.com/ChunHuangPhy/CompactObject.git
   cd CompactObject
   conda env create --file environment.yml
   conda activate CompactObject
   python -m pip install --no-deps --no-build-isolation -e .

The final command installs only the CompactObject source package. The
``--no-deps`` and ``--no-build-isolation`` options are intentional because the
runtime and build dependencies have already been installed by Conda. Using
``python -m pip`` also guarantees that pip belongs to the Python interpreter in
the active environment.

Do not use ``conda create -n CompactObject`` without listing Python. An empty
Conda environment has no local Python or pip executable, so a subsequent
``pip`` command can resolve to the system installation and fail with an
``externally-managed-environment`` error.

To update an existing environment after ``environment.yml`` changes, recreate
it:

.. code-block:: bash

   conda env remove --name CompactObject
   conda env create --file environment.yml

Dependency policy
-----------------

``pyproject.toml`` is the canonical dependency declaration for pip builds and
published package metadata. It lists direct, unpinned runtime dependencies so
pip can select a mutually compatible set. ``environment.yml`` is the
corresponding Conda-native environment and is tested separately in continuous
integration. Exact package versions belong in a generated lock file for a
specific reproducible analysis, not in the library's installation
instructions.

.. _optional-fastrmf-installation:

Optional FastRMF dependency
---------------------------

``CompactObject-TOV`` can use ``NumbaMinpack`` for the accelerated
``EOSgenerators.fastRMF_EoS`` path. This dependency is optional; the standard
RMF and DDH implementations work without it.

The `NumbaMinpack release on PyPI
<https://pypi.org/project/NumbaMinpack/>`_ (0.1.3, November 2021) can fail to
build with newer CMake versions. The `upstream GitHub source
<https://github.com/Nicholaswogan/NumbaMinpack>`_ contains updated CMake
requirements that are not included in that release. If the PyPI installation
fails with a CMake compatibility error, install from GitHub as described below,
following the recommendation in `CompactObject issue #137
<https://github.com/ChunHuangPhy/CompactObject/issues/137>`_.

NumbaMinpack requires a Fortran compiler. First install the compiler, CMake,
and Git. On macOS:

.. code-block:: bash

   brew install gcc cmake git

On Debian or Ubuntu:

.. code-block:: bash

   sudo apt-get update
   sudo apt-get install build-essential gfortran cmake git

After installing CompactObject, activate the same Python environment used for
CompactObject and the notebooks (whether managed by ``venv``, pyenv, or Conda).
Install NumbaMinpack directly from GitHub:

.. code-block:: bash

   python -m pip install --force-reinstall --no-deps "git+https://github.com/Nicholaswogan/NumbaMinpack.git"

Alternatively, clone the upstream repository and install the local source:

.. code-block:: bash

   git clone https://github.com/Nicholaswogan/NumbaMinpack.git
   python -m pip install --force-reinstall --no-deps ./NumbaMinpack

``--force-reinstall`` replaces any existing PyPI installation, since the GitHub
source currently has the same package version. ``--no-deps`` preserves the
NumPy and Numba versions already installed with CompactObject; pip still
installs the isolated build dependencies declared by NumbaMinpack.

Verify the installation in that environment:

.. code-block:: bash

   python -c "from NumbaMinpack import lmdif, hybrd, minpack_sig; print('NumbaMinpack import OK')"

Restart the notebook kernel after installation and ensure it uses this Python
environment.

From a CompactObject source checkout, the ``fast`` extra declares this optional
dependency but still resolves it from PyPI by default:

.. code-block:: bash

   python -m pip install -e ".[fast]"

Install the GitHub source first if you encounter the CMake error. A subsequent
installation of the ``fast`` extra will accept the installed NumbaMinpack. For
documentation and notebook development with the optional accelerated path:

.. code-block:: bash

   python -m pip install -e ".[docs,notebooks,fast]"

NumbaMinpack builds a native library, and the GitHub source may still fail on
some compiler, Python, or platform combinations. See, for example, upstream
`issue #10 <https://github.com/Nicholaswogan/NumbaMinpack/issues/10>`_. Report
remaining build or platform problems to the `upstream issue tracker
<https://github.com/Nicholaswogan/NumbaMinpack/issues>`_.

If the native build remains unavailable, omit the ``fast`` extra and use:

.. code-block:: bash

   python -m pip install -e ".[notebooks]"

The example notebook will skip only the FastRMF benchmark. The standard RMF
and DDH examples do not require NumbaMinpack.
