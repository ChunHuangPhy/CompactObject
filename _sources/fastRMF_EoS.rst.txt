.. _Fast_RMF_EOS:

RMF EOS solver (numba speed up version)
========================================

Functions to compute RMF Equation of state from given parameters. (numba speed up)

This module requires the optional ``NumbaMinpack`` dependency. Follow
:ref:`optional-fastrmf-installation` for compiler prerequisites, installation
from the upstream GitHub source when the PyPI build fails with a CMake error,
and an import check. Install it in the same Python environment as CompactObject
before importing this module.

From a CompactObject source checkout, the ``fast`` extra declares this
dependency but resolves it from PyPI unless it is already installed:

.. code-block:: bash

   python -m pip install -e ".[fast]"

The standard RMF implementation remains available from :ref:`RMF_EOS` without
this optional dependency.
