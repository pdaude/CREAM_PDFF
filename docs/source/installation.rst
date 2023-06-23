Installation instructions
=========================

We strongly recommend that you install CREAM_PDFF in a virtual environment !

.. code-block:: console

    git clone https://github.com/pdaude/CREAM_PDFF.git
    pip install -e CREAM_PDFF


Fatty_RIOT_GC repository
************************

#. Run *setup_berglund_QPBO.m* keeping it in the same folder as *test_berglund_QPBO.m*
    * Adds the parent folder of *test_berglund_QPBO.m* to the MATLAB path
    * Adds the parent folder of *fw_i3cm1i_3pluspoint_berglund_QPBO.m* to the MATLAB path
#. LINUX users need to unzip the file *./berglund/QPBO/DixonApp/LINUX/DixonApp_LINUX.exe.zip* (GITHUB does not allow files more than 100 MB in size)

.. note:: Information extracted from the original repository

VLGCA repository
****************

To convert CUDA file to mex file

.. code-block:: console

    cd fieldmap-graph-cut/code/gandalf/
    mexcuda residualcalculation_cuda.cu
