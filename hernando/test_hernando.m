clear all 
clc
close all
%Add Paths :
% with subfolders


% Algorithm-specific parameters
% algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
% algoParams.range_r2star = [0 100]; % Range of R2* values
% algoParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
% algoParams.range_fm = [-400 400]; % Range of field map values
% algoParams.NUM_FMS = 301; % Number of field map values to discretize
% algoParams.NUM_ITERS = 40; % Number of graph cut iterations
% algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
% algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
% algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
% algoParams.lambda = 0.05; % Regularization parameter
% algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
% algoParams.TRY_PERIODIC_RESIDUAL = 0;


algoParams=ReadYaml('../CREAM_PDFF/algoParams/Hernando_algoParams.yml');

[species ,FWspectrum]= setupModelParams(modelParams);
algoParams.gyro=FWspectrum.gyro;
algoParams.species =species;

img=load('../CREAM_PDFF/simu.mat');

[params, sse] = hernando_main(img.imDataParams,algoParams);
save('hernando_simu_SNR50.mat','params','sse'); 