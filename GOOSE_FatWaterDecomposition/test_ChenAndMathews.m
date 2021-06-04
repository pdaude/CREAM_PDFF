%% test_ChenAndMathews
%%
%% Test fat-water algorithms on (ISMRM Challege's) data
%%
%% Author: Pierre Daude

clear all 
clc
close all
%Add Paths :
%/home/pdaude/Projet_Python/CREAM_PDFF/Bydder with subfolders

%algoParams=ReadYaml('/home/pdaude/Projet_Python/CREAM_PDFF/algoParams/GOOSE_algoParams.yml');
img=load('/home/pdaude/Projet_Python/CREAM_PDFF/fw_i3cm1i_3pluspoint_berglund_QPBO/test_cases/01.mat');

% General parameters

% % General parameters
% algoParams.species(1).name = 'water';
% algoParams.species(1).frequency = 0;
% algoParams.species(1).relAmps = 1;
% algoParams.species(2).name = 'fat';
% algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
% algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
% 
% Numlayers = 100; % number of layers in the graph
% algoParams.gyro=42.57747892;
% gyro = algoParams.gyro%42.58;
% freqRange = gyro*img.imDataParams.FieldStrength*8; %search range is [-8ppm, 8ppm];
% t = linspace(-freqRange, freqRange, Numlayers);
% gridspacing = t(2)-t(1);
% 
% %% Set recon parameters
% % General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
algoParams.gyro=42.57747892;
% 
% % Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 30]; % Range of R2* values
algoParams.NUM_R2STARS = 16; % Numbre of R2* values for quantization
algoParams.NUM_FMS = 100; % Number of field map values to discretize also Num of Layers

algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 1; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;


[params,sse]=GOOSE_main(img.imDataParams,algoParams);
save('output_B0NICE.mat','params','sse'); 