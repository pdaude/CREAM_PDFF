clear all 
clc
close all
%Add Paths :
%/home/pdaude/Projet_Python/CREAM_PDFF/Bydder with subfolders

algoParams=ReadYaml('../CREAM_PDFF/algoParams/Bydder_algoParams.yml');

modelParams =ReadYaml('../CREAM_PDFF/modelParams/CustommodelParams.yml');

[species ,FWspectrum]= setupModelParams(modelParams);
algoParams.gyro=FWspectrum.gyro;
algoParams.species =species;

img=load('../CREAM_PDFF/simu.mat');

tic
[params, sse,sse1,YM] = Bydder_main(img.imDataParams,algoParams);
toc
save('output_Bydder_simu.mat','params','sse');
