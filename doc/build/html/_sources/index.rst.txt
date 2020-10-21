rpCofactors's Documentation
===========================

Indices and tables
##################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Introduction
############

.. _RetroRules: https://retrorules.org/
.. _RetroPath2.0: https://github.com/Galaxy-SynBioCAD/RetroPath2
.. _rp2paths: https://github.com/Galaxy-SynBioCAD/rp2paths
.. _rpSBML: https://github.com/Galaxy-SynBioCAD/rpBase
.. _rpBase: https://github.com/Galaxy-SynBioCAD/rpBase
.. _rpCache: https://github.com/Galaxy-SynBioCAD/rpCache
.. _rpReader: https://github.com/Galaxy-SynBioCAD/rpReader
.. _rpUnicity: https://github.com/Galaxy-SynBioCAD/rpUnicity

Welcome rpCofactors's documentation. This tool adds cofactors to mono-component reaction rules from RetroRules_ and must be used after rpReader_.

.. note::
   This projects relies on rpUnicity_ project to filter any duplicates that may occur due to differen reaction rules having the same cofactors.

To build the docker you must build a rpBase_ and rpCache_ docker, and then you can use the following command:

.. code-block:: bash

   docker build -t brsynth/rpcofactors-standalone:v2 .

You can run the docker using the following command to parse the rpReader_:

.. code-block:: bash

   python run.py -input test/test_rpReader.tar -output test/test_rpCofactors.tar -input_format tar

API
###

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. currentmodule:: rpTool

.. autoclass:: rpCofactors
    :show-inheritance:
    :members:
    :inherited-members:

.. currentmodule:: rpToolServe

.. autoclass:: main
    :show-inheritance:
    :members:
    :inherited-members:

.. autoclass:: main_extrules
    :show-inheritance:
    :members:
    :inherited-members:

.. autoclass:: runSingleSBML
    :show-inheritance:
    :members:
    :inherited-members:

.. autoclass:: runCofactors_mem
    :show-inheritance:
    :members:
    :inherited-members:

.. autoclass:: runCofactors_hdd
    :show-inheritance:
    :members:
    :inherited-members:

.. currentmodule:: run

.. autoclass:: main
    :show-inheritance:
    :members:
    :inherited-members:

