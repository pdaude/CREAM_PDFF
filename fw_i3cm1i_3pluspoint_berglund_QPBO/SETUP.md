# fw_i3cm1i_3pluspoint_berglund_QPBO

Requirements:
=============
* Mac OS X (tested with Lion) or Windows (tested with 32-bit WinXP) or i386 LINUX (tested on Ubuntu 8.04 LTS)
* Tested with MATLAB R0210a, MATLAB 2012a and MATLAB 2013a

Installation:
=============
* Run `setup_berglund_QPBO.m` keeping it in the same folder as `test_berglund_QPBO.m`
  - Adds the parent folder of `test_berglund_QPBO.m` to the MATLAB path
  - Adds the parent folder of `fw_i3cm1i_3pluspoint_berglund_QPBO.m` to the MATLAB path
* LINUX users need to unzip the file `./berglund/QPBO/DixonApp/LINUX/DixonApp_LINUX.exe.zip` (GITHUB does not allow files more than 100 MB in size)

Properties:
======
- Image-space
- 2 species (water-fat)
- Complex-fitting
- Multi-peak fat (pre-calibrated)
- R2* (given >3 echoes)
- Independent water/fat phase
- Requires 3 or more uniformly spaced echoes

Usage:
======
`outParams = fw_i3cm1i_3pluspoint_berglund_QPBO( imDataParams, algoParams );`

### Basic Input: ISMRM fat-water toolbox structure

`imDataParams`

- `imDataParams.images `               : acquired images, array of size [nx, ny, nz, ncoils, nTE]
- `imDataParams.TE`                    : echo times (in seconds), vector of length nTE
- `imDataParams.FieldStrength`         : (in Tesla), (default 1.5)
- `imDataParams.PrecessionIsClockwise` : ==1 is clockwise, ~=1 is counterclockwise, (default 1) 
- `imDataParams.mask`                  : logical [nx, ny, nz], ( default true([nx, ny, nz]) )

`algoParams`

- `algoParams.species(ii).name`        : name of species ii (string)
- `algoParams.species(ii).frequency`   : frequency shift in ppm of each peak within species ii
- `algoParams.species(ii).relAmps`     : relative amplitude of each peak within species ii

### Advanced Input: options specific to QPBO algorithm

- `algoParams.decoupled_estimation = true;`          : flag for decoupled R2 estimation
- `algoParams.Fibonacci_search = true;`              : flag for Fibonacci search
- `algoParams.B0_smooth_in_stack_direction = false;` : flag for B0 smooth in stack direction
- `algoParams.multigrid = true;`                     : flag for multi-level resolution pyramid
- `algoParams.estimate_R2 = true;`                   : flag to estimate R2star
- `algoParams.verbose = false;`                      : flag for verbose status messages (default false)
- `algoParams.process_in_3D = true;`                 : flag to process in 3D (default true)
- `algoParams.R2star_calibration = false;`           : flag to perform R2* calibration (default false)
- `algoParams.ICM_iterations = 2;`                   : ICM iterations
- `algoParams.num_B0_labels = 100;`                  : number of discretized B0 values
- `algoParams.mu = 10;`                              : regularization parameter
- `algoParams.R2_stepsize = 1;`                      : R2 stepsize in s^-1
- `algoParams.max_R2 = 120;`                         : maximum R2 in s^-1
- `algoParams.max_label_change = 0.1;`               : 
- `algoParams.fine_R2_stepsize = 1.0;`               : Fine stepsize for discretization of R2(*) [s^-1] (used in decoupled fine-tuning step)
- `algoParams.coarse_R2_stepsize = 10.0;`            : Coarse stepsize for discretization of R2(*) [s^-1] (used for joint estimation step, value 0.0 --> Decoupled estimation only
- `algoParams.water_R2 = 0.0;`                       : Water R2 [sec-1]
- `algoParams.fat_R2s = zeros(1,9);`                 : fat peak R2s [s^-1]
- `algoParams.R2star_calibration_max = 800;`         : max R2* to use for calibration [s^-1] (default 800)
- `algoParams.R2star_calibration_cdf_threshold;`     : threshold for R2* calibration cumulative density function [0,1] (default 0.98)

### Output:
- `outParams.species(ii).name`       : name of the species (taken from algoParams)
- `outParams.species(ii).amps`       : estimated water/fat images, size [nx,ny,nz] 
- `outParams.fieldmap`               : field map (in Hz, size [nx,ny,nz]), fieldmap is NOT unwrapped
- `outParams.r2starmap`              : R2* map (in s^{-1}, size [nx,ny,nz])
- `outParams.status`                 : QPBO status map (size [nx,ny,nz])
- `outParams.residual`               : residual error map (size [nx,ny,nz])
- `outParams.mask`                   : mask (size [nx,ny,nz])
- `outParams.execution_time_seconds` : execution time in seconds

### Example Species Description
- `algoParams.species(1).name = 'water';`
- `algoParams.species(1).frequency = 4.70;`
- `algoParams.species(1).relAmps = 1;`
- `algoParams.species(2).name = 'fat';`
- `algoParams.species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29];` : 9-peak model
- `algoParams.species(2).relAmps   = [  88,  642,   58,   62,   58,    6,   39,   10,   37];` : Hamilton G, et al. NMR Biomed. 24(7):784-90, 2011. PMID: 21834002
  
Example:
========
* `test_berglund_QPBO.m`
  - Tests `fw_i3cm1i_3pluspoint_berglund_QPBO.m` on all 10 Phase I and 7 Phase II cases from the [2012 ISMRM Challenge](http://www.ismrm.org/challenge/node/18)
