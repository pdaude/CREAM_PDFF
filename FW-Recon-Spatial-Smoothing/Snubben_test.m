clear all 
clc
close all
%Add Paths :
%/home/pdaude/Projet_Python/CREAM_PDFF/Bydder with subfolders
%addpath(['.' filesep 'qpboMex-master'])

algoParams=ReadYaml('/home/pdaude/Projet_Python/CREAM_PDFF/algoParams/Snubben_algoParams.yml');
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 4.7;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [0.9,1.3,2.1,2.76,4.31,5.3];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
algoParams.gyro=42.57747892;

param.fatCS = [5.3,4.31,2.76,2.1,1.3,0.9];%The chemical shifts of the fat peaks, Unit: ppm
param.relAmps = [0.048,0.039,0.004,0.128,0.693,0.087];%The relative weights of the fat peaks
param.watCS = 4.7;%The chemical shift of water, Unit: ppm



img=load('/home/pdaude/Projet_Python/CREAM_PDFF/fw_i3cm1i_3pluspoint_berglund_QPBO/test_cases/01.mat');

[params, sse] = Snubben_main(img.imDataParams,algoParams);
save('/home/pdaude/Projet_Python/CREAM_PDFF/Bydder/output_Bydder_01.mat','params','sse');
