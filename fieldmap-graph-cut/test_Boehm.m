%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,nz,1,nTE]
%%   - imDataParams.TE: echo times (in seconds)
%%   - imDataParams.FieldStrength: (in Tesla)
%%
%%   - algoParams.species(ii).name = name of species ii (string)
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%%   Example
%%      - algoParams.species(1).name = 'water' % Water
%%      - algoParams.species(1).frequency = [0] 
%%      - algoParams.species(1).relAmps = [1]   
%%      - algoParams.species(2).name = 'fat' % Fat
%%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%% 
%%   - algoParams.range_r2star = [0 0]; % Range of R2* values
%%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
%%   - algoParams.range_fm = [-400 400]; % Range of field-map values
%%   - algoParams.NUM_FMS = 301; % Number of field-map values to discretize
%%
%% Output: structure waterfatparams 
%%   - waterfatparams.water: water image
%%   - waterfatparams.fat: fat image
%%   - waterfatparams.: fat image
%%   - waterfatparams.fat: fat image
%%   - outParams.r2starmap: R2* map (in s^{-1}, size [nx,ny,nz])
%%   - outParams.fieldmap: field-map (in Hz, size [nx,ny,nz])
%%
%% Author: Pierre Daude

clear all 
clc
close all
%Add Paths :
% with subfolders

algoParams=ReadYaml('../CREAM_PDFF/algoParams/Boehm_algoParams.yml');

modelParams =ReadYaml('../CREAM_PDFF/modelParams/CustommodelParams.yml');
[species ,FWspectrum]= setupModelParams(modelParams);
algoParams.gyro=FWspectrum.gyro;
algoParams.species =species;

img=load('../CREAM_PDFF/simu.mat');

[params, sse] = GC_main(img.imDataParams,algoParams);
save('output_Boehm_simu.mat','params','sse'); 