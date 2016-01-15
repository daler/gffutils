.. _api:


.. rst-class:: html-toggle

API
===

`gffutils`
----------

.. automodule:: gffutils
    :members:

Create a database
-----------------


.. currentmodule:: gffutils.create

.. autosummary::
    :toctree: autodocs
    :nosignatures:

    create_db

Interact with a database
------------------------
First, connect to an existing database:

.. currentmodule:: gffutils.interface

.. autosummary::
    :toctree: autodocs
    :nosignatures:

    FeatureDB

Then, use the methods of :class:`FeatureDB` to interact:


.. autosummary::
    :toctree: autodocs
    :nosignatures:

    FeatureDB.children
    FeatureDB.parents
    FeatureDB.schema
    FeatureDB.features_of_type
    FeatureDB.count_features_of_type
    FeatureDB.all_features
    FeatureDB.execute
    FeatureDB.featuretypes
    FeatureDB.region
    FeatureDB.iter_by_parent_childs

Modify a :class:`FeatureDB`:

.. autosummary::
    :toctree: autodocs
    :nosignatures:

    FeatureDB.update
    FeatureDB.delete
    FeatureDB.add_relation
    FeatureDB.set_pragmas

Operate on features:

.. autosummary::
    :toctree: autodocs
    :nosignatures:

    FeatureDB.interfeatures
    FeatureDB.children_bp
    FeatureDB.merge
    FeatureDB.create_introns
    FeatureDB.bed12




Feature objects
---------------
Most :class:`FeatureDB` methods return :class:`Feature` objects:


.. currentmodule:: gffutils.feature

.. autosummary::
    :toctree: autodocs
    :nosignatures:

    Feature

You can extract the sequence for a feature:

.. autosummary::
    :toctree: autodocs
    :nosignatures:

    Feature.sequence

Creating a :class:`Feature` object:

.. autosummary::
    :toctree: autodocs
    :nosignatures:

    feature_from_line

Integration with other tools
----------------------------
.. currentmodule:: gffutils

.. autosummary::
    :toctree: autodocs
    :nosignatures:

    gffutils.biopython_integration.to_seqfeature
    gffutils.biopython_integration.from_seqfeature
    gffutils.pybedtools_integration.tsses
    gffutils.pybedtools_integration.to_bedtool



Utilities
---------

.. autosummary::
    :toctree: autodocs
    :nosignatures:

    gffutils.helpers.asinterval
    gffutils.helpers.merge_attributes
    gffutils.helpers.sanitize_gff_db
    gffutils.helpers.annotate_gff_db
    gffutils.helpers.infer_dialect
    gffutils.helpers.example_filename
    gffutils.inspect.inspect
