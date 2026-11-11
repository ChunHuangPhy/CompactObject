Releasing
=========

PyPI releases are immutable. Changes to package metadata, including optional
extras such as ``test``, require a new version and cannot be added to an
existing release.

Trusted Publisher setup
-----------------------

CompactObject uses the dedicated ``.github/workflows/publish.yaml`` workflow
for PyPI publishing. A PyPI project owner must authorize it once:

1. Open the ``CompactObject-TOV`` project on PyPI.
2. Select **Manage**, then **Publishing**.
3. Add a GitHub Actions trusted publisher with these values:

   - Owner: ``ChunHuangPhy``
   - Repository: ``CompactObject``
   - Workflow name: ``publish.yaml``
   - Environment name: ``pypi``

The environment name must match exactly. Configure the ``pypi`` environment in
the GitHub repository with required maintainer approval before deployment.
PyPI's `Trusted Publisher documentation
<https://docs.pypi.org/trusted-publishers/adding-a-publisher/>`_ describes the
one-time project configuration.

This setup uses a short-lived OpenID Connect credential. It does not require a
PyPI API token to be stored in the repository.

Publishing a version
--------------------

Before creating a release:

1. Update the version in ``pyproject.toml`` and ``docs/source/conf.py``.
2. Run ``python -m pytest``.
3. Build and validate the distributions:

   .. code-block:: bash

      python -m build
      python -m twine check dist/*

4. Commit and push the release changes.
5. Create a Git tag whose name is ``v`` followed by the package version.
6. Publish a GitHub release from that tag.

For version 2.1.1, the expected tag is ``v2.1.1``. Publishing that GitHub
release starts the PyPI workflow. The workflow verifies that the tag and
``pyproject.toml`` version match, builds and checks both distributions, and
uploads them to PyPI.
