clear all 
clc
close all

% param.fatCS = [5.3,4.31,2.76,2.1,1.3,0.9];%The chemical shifts of the fat peaks, Unit: ppm
% param.relAmps = [0.048,0.039,0.004,0.128,0.693,0.087];%The relative weights of the fat peaks
% param.watCS = 4.7;%The chemical shift of water, Unit: ppm

%Add Paths :
%addpath(['.' filesep 'qpboMex-master'])

algoParams=ReadYaml('../CREAM_PDFF/algoParams/Snubben_algoParams.yml');
modelParams =ReadYaml('../CREAM_PDFF/modelParams/CustommodelParams.yml');
[species ,FWspectrum]= setupModelParams(modelParams);
algoParams.gyro=FWspectrum.gyro;
algoParams.species =species;

img=load('../CREAM_PDFF/simu.mat');

[params, sse] = Snubben_main(img.imDataParams,algoParams);
save('output_snubben_simu.mat','params','sse');
