*Copyright (c) 2016–2020 Johan Berglund*

*FWQPBO is distributed under the terms of the GNU General Public License*

*This program is free software: you can redistribute it and/or modify*
*it under the terms of the GNU General Public License as published by*
*the Free Software Foundation, either version 3 of the License, or*
*(at your option) any later version.*

*This program is distributed in the hope that it will be useful,*
*but WITHOUT ANY WARRANTY; without even the implied warranty of*
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the*
*GNU General Public License for more details.*

*You should have received a copy of the GNU General Public License*
*along with this program.  If not, see <http://www.gnu.org/licenses/>.*

ABOUT
-------------------------------------------------------------------------------
FWQPBO is a command-line tool for MRI [chemical shift based fat-water 
separation with B0-correction based on QPBO graph cuts
](https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.26479). 
Input paramaters are provided by human readable configuration files. 
See example configuration files provided with FWQPBO for details. 
Input data can be in DICOM format, or [MATLAB format according to the 
ISMRM 2012 challenge](http://challenge.ismrm.org/node/14).
FWQPBO is written in Python.

HOW TO USE
-------------------------------------------------------------------------------
First install required packages, see dependencies.
To use the command-line tool, type 
```
./main.py -h.
```
To use as a Python script, see example file [./demo.py](demo.py). The demo
reconstructs [data from the ISMRM 2012 challenge (data download required)
](http://challenge.ismrm.org/node/4).

HOW TO CITE
-------------------------------------------------------------------------------
Berglund J and Skorpil M. *Multi-scale graph-cut algorithm for efficient water-
fat separation*. Magn Reson Med, 78(3):941-949, 2017. [doi: 10.1002/mrm.26479]

DEPENDENCIES
-------------------------------------------------------------------------------
See [./environment.yml](environment.yml)

CONTACT INFORMATION
-------------------------------------------------------------------------------
Johan Berglund, Ph.D.  
Dept. of Clinical Neuroscience  
Karolinska Institutet,  
Stockholm, Sweden  
johan.berglund@ki.se  
