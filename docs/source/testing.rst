Testing
*******

The automated tests live in the repository's ``tests/`` directory and are not
included in the installed wheel. Clone the source repository before running
them.

Fresh setup with venv
---------------------

Create a source checkout and isolated Python environment:

.. code-block:: bash

   git clone https://github.com/ChunHuangPhy/CompactObject.git
   cd CompactObject
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install -e ".[test]"

On Windows, activate the environment with ``.venv\Scripts\activate``.

Run the complete suite from the repository root:

.. code-block:: bash

   python -m pytest

The ``test`` extra was added in CompactObject 2.1.1. It installs pytest in
addition to the package's runtime dependencies.

Conda setup
-----------

The supplied Conda environment already includes pytest:

.. code-block:: bash

   git clone https://github.com/ChunHuangPhy/CompactObject.git
   cd CompactObject
   conda env create --file environment.yml
   conda activate CompactObject
   python -m pip install --no-deps --no-build-isolation -e .
   python -m pytest

Running selected tests
----------------------

Pass a file path to run one test module:

.. code-block:: bash

   python -m pytest tests/test_ddh_char23.py

Pass a pytest node ID to run one test:

.. code-block:: bash

   python -m pytest tests/test_ddh_char23.py::DDHChar23Tests::test_builtin_char23_matches_user_defined_form

Use ``-vv`` for more detailed output:

.. code-block:: bash

   python -m pytest -vv

Troubleshooting the test extra
------------------------------

If pip reports that ``compactobject-tov 2.1 does not provide the extra
'test'``, the environment or source checkout is still using the older 2.1
metadata. Update the source checkout and reinstall it:

.. code-block:: bash

   git switch main
   git pull --ff-only
   python -m pip install -e ".[test]"

Pytest can also be installed explicitly while working from a source checkout:

.. code-block:: bash

   python -m pip install -e .
   python -m pip install pytest
   python -m pytest

Confirm which CompactObject release is installed with:

.. code-block:: bash

   python -c "from importlib.metadata import version; print(version('CompactObject-TOV'))"

Continuous integration
----------------------

The pytest configuration is stored in ``pyproject.toml``. Runtime warnings are
treated as failures so numerical domain errors cannot pass unnoticed. GitHub
Actions runs the complete suite automatically for pushes and pull requests.

Polytrope tutorial checks
-------------------------

``test_Inference_polytrope.ipynb`` separates the inexpensive EOS unit check
from stellar integrations and full inference. Its default unit check draws
three EOS curves without calling the TOV solver. The optional mass-radius
survey is enabled with ``RUN_MASS_RADIUS_CHECK=True`` and prints progress for
each curve.

The default UltraNest interface check uses eight live points and a coarse
six-point central-density grid over log10 density 14.6–15.6 with a small
call/iteration budget. This restricted demonstration prior checks execution
only; its output is not a converged posterior or evidence estimate.
The full inference is enabled separately with ``RUN_FULL_INFERENCE=True`` and
uses the 50-point density grid over log10 density 14.3–15.6. Individual TOV integrations can remain slow,
so a call budget is not a wall-clock timeout.

The documentation workflow executes the default notebook checks with
``jupyter nbconvert`` and excludes cells tagged ``long-running``. Regression
tests exercise both the tutorial and documentation notebook copies, checking
EOS units, the density grid used by the prior, empty stellar sequences, and
the absence of TOV calls in the default unit check:

.. code-block:: bash

   python -m pytest tests/test_polytrope_notebook.py
