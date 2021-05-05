param.nR2 = 145;%This is the number of different values of R2* that will be considered
param.R2max = 144.0;%This is the maximal value of R2* that will be considered, Unit: s^-1
param.R2cand = 40.;%This is the initial guess of R2*, Unit: s^-1
param.nB0 = 100;%This is the number of off-resonance frequencies that will be considered
param.use3D = 1;%set to 0 for 2D, 1 for 3D
param.nICMiter = 10;%Number of ICM iterations
param.maxICMupdate = round(param.nB0/10);%Restriction of the range of the values ICM will consider in each iteration
param.sigma = 2^.75;%Value of sigma, Unit: mm^-1
param.lambda = 10;%Initial value of lambda

param.fatCS = [5.3,4.31,2.76,2.1,1.3,0.9];%The chemical shifts of the fat peaks, Unit: ppm
param.relAmps = [0.048,0.039,0.004,0.128,0.693,0.087];%The relative weights of the fat peaks
param.watCS = 4.7;%The chemical shift of water, Unit: ppm
param.gyro = 42.576;%Gyromagnetic ratio, Unit: MHz

voxelsize = [1.5 1.5 5];%Voxel sizes, [dx dy dz], Unit: mm

load_file = '01.mat';%Should be in the same format as in the 2012 ISMRM Challenge is on water-fat reconstruction: http://challenge.ismrm.org/node/14
tic
[water, fat, FF, R2map, B0map, label, lambda_map] = doReconGauss( load_file, voxelsize, param );
toc