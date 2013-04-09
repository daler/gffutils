.. _api:

.. currentmodule:: gffutils

.. rst-class:: html-toggle

API
===
Create a database
-----------------
.. autosummary::
    :toctree: autodocs

    gffutils.create_db

Interact with a database
------------------------
First, connect to an existing database:

.. autosummary::
    :toctree: autodocs

    gffutils.FeatureDB

Then, use the methods of :class:`FeatureDB` to interact:

.. autosummary::
    :toctree: autodocs

    gffutils.FeatureDB.children
    gffutils.FeatureDB.parents
    gffutils.FeatureDB.schema
    gffutils.FeatureDB.features_of_type
    gffutils.FeatureDB.count_features_of_type
    gffutils.FeatureDB.all_features
    gffutils.FeatureDB.execute
    gffutils.FeatureDB.featuretypes
    gffutils.FeatureDB.region


Many :class:`FeatureDB` methods return :class:`Feature` objects:

.. autosummary::
    :toctree: autodocs

    gffutils.Feature


If you'd like to parse GFF/GTF files yourself, try the :class:`Parser` class:

.. autosummary::
    :toctree: autodocs

    gffutils.Parser
    gffutils.Parser.show_warnings
