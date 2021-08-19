%% setup algoParams
clc
close all
clear all

algoParams=ReadYaml('../CREAM_PDFF/algoParams/Berglund_algoParams.yml');
modelParams =ReadYaml('../CREAM_PDFF/modelParams/CustommodelParams.yml');
[species ,FWspectrum]= setupModelParams(modelParams);
algoParams.gyro=FWspectrum.gyro;
algoParams.species =species;

img=load('../CREAM_PDFF/simu.mat');
[params,sse,outParams]=Berglund_main(img.imDataParams,algoParams);
save('output_Berglund_simu.mat','params','sse'); 