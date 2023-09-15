How does it work ?
==================

.. image:: howdoesitwork.JPG



Each fat-water separation algorithm have three input structures (imDataParams modelParams and algoParams) and one output structure (outParams). 

Input Data (imDataParams)
*************************
The input data should be inside MATLAB file (.mat ) and contain the following fields:

- imDataParams.images : the data, array of dimensions nx X ny X nz X ncoils X nte
- imDataParams.FieldStrength : field strength (in Tesla)
- imDataParams.TE : echo times (in s) ; vector of length nte
- imDataParams.PrecessionIsClockwise (1/-1)
- imDataParams.voxelSize : size of the voxel of the image vector of length 3

Algorithm parameters (algoParams)
*********************************
The algorithm parameters are loaded from a YAML file (.yml). In the algoParams folder, you can find an example of YAML file for each fat-water separation algorithm tested in this toolbox. 

Input fat spectrum model (modelParams)
**************************************
The input fat spectrum model are also loaded from a YAML file (.yml). In the modelParams folder, each YAML file corresponds to a specific fat spectrum model. 
