clear all 
clc
close all
%Add Paths :
%/home/pdaude/Projet_Python/CREAM_PDFF/Bydder with subfolders

algoParams=ReadYaml('../CREAM_PDFF/algoParams/Bydder_algoParams.yml');

modelParams =ReadYaml('../CREAM_PDFF/modelParams/Hodson2008modelParams.yml');

[species ,FWspectrum]= setupModelParams(modelParams);
algoParams.gyro=FWspectrum.gyro;
algoParams.species =species;

img=load('../CREAM_PDFF/simulation/simuImDataParams/simuImDataParams_SNR100_TE9ID_Hodson2008.mat');

tic
[params, sse,sse1,YM] = Bydder_main(img.imDataParams,algoParams);
toc
save('output_Bydder_simu.mat','params','sse','sse1');
