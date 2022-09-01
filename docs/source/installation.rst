
Installation
============

scmorph is currently available through `Github <https://github.com/jeskowagner/scmorph>`_ and can be installed with ``pip`` (or any equivalent):

.. code-block:: bash

   pip install git+https://github.com/jeskowagner/scmorph/

Alternatively, you may also install the dependencies via ``conda``. Note that at the time of writing, the package itself is not yet released on conda and will therefore be installed through ``pip``.

.. code-block:: bash

   tmpfile=$(mktemp)

   # Download dependencies list
   wget -O $tmpfile https://raw.githubusercontent.com/jeskowagner/scmorph/8cd6111e0cf0ad26a29b7ad663cf6dcc312128f7/env.yml?token=GHSAT0AAAAAABWL2SDH24EA54CW2NLO56JIYXNJQMQ

   # Create environment from depencies
   conda env create --file $tmpfile -n scmorph
   conda activate scmorph

   # Install scmorph
   pip install git+https://github.com/jeskowagner/scmorph/

In the future we will release scmorph on conda and PyPI, so watch this space for updates!

.. raw:: html

   <!--
   The package is available through [PyPI](https://pypi.org/project/scmorph/) and can be installed with `pip` (or any equivalent):
   ```bash
   pip install scmorph
   ```
   -->
