%% test_hernando_111110
%%
%% Test fat-water algorithms on (Peter Kellman's) data
%%
%% Author: Diego Hernando
%% Date created: August 18, 2011
%% Date last modified: February 29, 2011

% Add to matlab path
BASEPATH = './';
addpath([BASEPATH 'common/']);
addpath([BASEPATH 'graphcut/']);
addpath([BASEPATH 'descent/']);
addpath([BASEPATH 'mixed_fitting/']);
addpath([BASEPATH 'create_synthetic/']);
addpath([BASEPATH 'matlab_bgl/']);



%% Load some data
foldername = [BASEPATH '../../fwtoolbox_v1_data/kellman_data/'];
fn = dir([foldername '*.mat']);
file_index = ceil(rand*length(fn));
disp([foldername fn(file_index).name]);
load([foldername fn(file_index).name]);
imDataParams = data;
imDataParams.images = double(data.images);
% $$$ imDataParams.FieldStrength = 1.5;
% $$$ imDataParams.PrecessionIsClockwise = -1;

%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 100]; % Range of R2* values
algoParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
algoParams.range_fm = [-400 400]; % Range of field map values
algoParams.NUM_FMS = 301; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
THRESHOLD = 0.01;

%% Recon -- graph cut 
%% (Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large 
%% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.)
tic
  outParams = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, algoParams );
toc

DO_MIXED_FIT = 0;
if DO_MIXED_FIT > 0 
  try
    %% Recon -- mixed fit for phase error correction
    % Initialize mixed fitting to graph cut solution
    algoParams.fieldmap = outParams.fieldmap;
    algoParams.r2starmap = outParams.r2starmap;
    algoParams.NUM_MAGN = 1;
    algoParams.THRESHOLD = 0.04;
    algoParams.range_r2star = [0 200];

    % Do mixed fitting
    %% (Hernando D, Hines CDG, Yu H, Reeder SB. Addressing phase errors in fat-water imaging 
    %% using a mixed magnitude/complex fitting method. Magn Reson Med; 2011.)
    outParamsMixed = fw_i2xm1c_3pluspoint_hernando_mixedfit( imDataParams, algoParams )
  end
end


