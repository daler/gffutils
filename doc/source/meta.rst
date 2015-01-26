Meta-docs (docs about the docs)
===============================
These docs are generated using `sphinx <http://sphinx-doc.org/>`_ 

All code in the documentation doubles as doctests, so it is guaranteed to run
correctly.  In addition, since the API docs (:doc:`api`) include docstrings
from the source code itself, these tests are run as well when the documentation
is generated.

Building the docs
-----------------
* Install the requirements from the `dev-requirements.txt` file in the top dir
  of the repo::

    $ pip install -r dev-requirements.txt

* Navigate to the :file:`doc` folder in the source.

* Run doctests::

    $ make doctest

* Build HTML.  The results are in :file:`doc/build/html/`::

    $ make html
