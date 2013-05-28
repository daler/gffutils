Meta-docs (docs about the docs)
===============================
These docs are generated using `sphinx <http://sphinx-doc.org/>`_ and the
fantastic `Cloud Sphinx Theme <http://pythonhosted.org/cloud_sptheme/>`_.  All
code in the documentation doubles as doctests, so it is guaranteed to run
correctly.  In addition, since the API docs (:doc:`api`) include docstrings
from the source code itself, these tests are run as well when the documentation
is generated.

Building the docs
-----------------
* Get the requirements::

    $ pip install -U sphinx cloud_sptheme

* Navigate to the :file:`doc` folder in the source.

* Run doctests::

    $ make doctest

* Build HTML.  The results are in :file:`doc/build/html/`::

    $ make html

