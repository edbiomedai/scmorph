
Installation
============

scmorph is currently available through `PyPi <https://pypi.org/>`_ and can be installed with ``pip``:

.. code-block:: bash

   pip install scmorph

Alternatively, you may also install the dependencies via ``conda``.
Note that at the time of writing, the package itself is not yet released on conda and will therefore be installed through ``pip``.
On Unix system, you can install the dependencies with:

.. code-block:: bash

   tmpfile=$(mktemp)

   # Download dependencies list
   wget -O $tmpfile https://raw.githubusercontent.com/edbiomedai/scmorph/main/env.yml?token=GHSAT0AAAAAABWL2SDHOADDKGA65BGODCPKYYR4DXA

   # Create environment from depencies
   conda env create --file $tmpfile -n scmorph
   conda activate scmorph

   # Install scmorph
   pip install scmorph

In the future we will release scmorph on conda-forge, so watch this space for updates!
