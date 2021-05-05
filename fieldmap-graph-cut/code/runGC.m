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
%% Author: Christof Boehm
%% Date created: September 16, 2020
%% MATLAB R2017b

addpath(genpath('.'))
path2data = '../data/exampleData.mat';
load(path2data)

algoParams.species(1).name = 'water'; % Water
algoParams.species(1).frequency = [0] ;
algoParams.species(1).relAmps = [1]; 
algoParams.species(2).name = 'fat'; % Fat
algoParams.species(2).frequency = [-3.8000 -3.4000 -3.1000 -2.6800 -2.4600 -1.9500 -0.5000 0.4900 0.5900];
algoParams.species(2).relAmps = [0.0899 0.5834 0.0599 0.0849 0.0599 0.0150 0.0400 0.0100 0.0569];

algoParams.range_r2star = [0 500]; % Range of R2* values
algoParams.NUM_R2STARS = 26; % Numbre of R2* values for quantization
algoParams.range_fm = [-450 450]; % Range of field map values
algoParams.NUM_FMS = 401; % Number of field map values to discretize

%% check validaty of parameters
algoParams = checkParamsAndSetDefaults_GANDALF(imDataParams, algoParams);

%% calculate VARPRO residual and extract minima 
VARPROparams = calculateMinimaDirect(imDataParams, algoParams)

%% perform core graph cut function that outputs field-map, water, fat and r2*-map
waterFatParams = GANDALF(imDataParams, algoParams, VARPROparams)

%% plot
iz = ceil(size(waterFatParams.water, 3) / 2);
close all;
h = figure
ax = subplot(2, 2, 1)
colormap(ax, gray)
imagesc(abs(waterFatParams.water(:, :, iz)))
axis image
xlabel('water [a.u.]')
set(gca,'xtick',[])
set(gca,'ytick',[])

ax = subplot(2, 2, 2)
imagesc(abs(waterFatParams.fat(:, :, iz)))
colormap(ax, gray)
axis image
xlabel('fat [a.u.]')
set(gca,'xtick',[])
set(gca,'ytick',[])

ax = subplot(2, 2, 3)
imagesc(waterFatParams.fieldmap(:, :, iz))
colormap(ax, 'plasma')
axis image
caxis([-200, 200])
colorbar
xlabel('field-map [Hz]')
set(gca,'xtick',[])
set(gca,'ytick',[])

ax = subplot(2, 2, 4)
imagesc(waterFatParams.r2starmap(:, :, iz))
colormap(ax, 'magma')
axis image
colorbar
xlabel('R_2^* [Hz]')
set(gca,'xtick',[])
set(gca,'ytick',[])
