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

Optional FastRMF dependency
---------------------------

``CompactObject-TOV`` can use ``NumbaMinpack`` for the accelerated
``EOSgenerators.fastRMF_EoS`` path. This dependency is optional; the standard
RMF and DDH implementations work without it.

``NumbaMinpack`` requires a Fortran compiler. On macOS, install GCC before the
package:

.. code-block:: bash

   brew install gcc
   python -m pip install NumbaMinpack

On Debian or Ubuntu, install ``gfortran`` and CMake first:

.. code-block:: bash

   sudo apt-get install gfortran cmake
   python -m pip install NumbaMinpack

``NumbaMinpack`` builds a native library during installation. Some current
CMake, compiler, Python, and platform combinations still fail in that upstream
build even after the documented prerequisites are installed. For example,
upstream `issue #10
<https://github.com/Nicholaswogan/NumbaMinpack/issues/10>`_ reports a CMake
compatibility failure on a newer Linux/Python environment.

This does not block the standard CompactObject solvers. If the native build
fails, omit the ``fast`` extra and continue with:

.. code-block:: bash

   python -m pip install -e ".[docs,notebooks]"

The example notebook will skip only the FastRMF benchmark. Build and platform
problems with ``NumbaMinpack`` should be reported to its `upstream issue
tracker <https://github.com/Nicholaswogan/NumbaMinpack/issues>`_.

From a source checkout, the equivalent optional extra is:

.. code-block:: bash

   python -m pip install -e ".[fast]"

For documentation and notebook development with the optional accelerated path:

.. code-block:: bash

   python -m pip install -e ".[docs,notebooks,fast]"
