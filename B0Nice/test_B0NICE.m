clear all 
clc
close all
%Add Paths :
%/home/pdaude/Projet_Python/CREAM_PDFF/B0Nice/mfiles/mfile_B0NICEmpt
%/home/pdaude/Projet_Python/CREAM_PDFF/B0Nice/mfiles/PUROR2D_MatlabCode


% % algoParams.debug = 2; % 0 off all the figure; 1 on all; 2 only on the final
% % 
% % %set up default algoParams
% % %for initial B0 mapping
% % algoParams.FILTsize = 0;
% % algoParams.th4unwrap = 0; %prctile slice by slice
% % algoParams.th4supp = 20;%prctile slice by slice
% % algoParams.th4STAS = 40;%prctile slice by slice
% % algoParams.unwrap_mode = 2; % 2 for 2D PUROR; 3 for 3D PUROR
% % % for magnitude fitting
% % algoParams.FF_trl = [0.0 0.15 0.25 0.5 0.75 0.85 1.0]; 
% % %algoParams.FF_trl = 0:0.1:1.0;
% % algoParams.R2star_pt = 0:1:300; % unit: Hz
% % %algoParams.FFth4magWater = 0.3;
% % %---------------------------------------------
% % % phase error correction
% % algoParams.flag_debug = 0;% 0 off all the figure; 1 on all; 2 only on the final
% % algoParams.th4RegionDivision = pi/2;
% % algoParams.minAreaUsingPhaseGradient = 20; % 
% % algoParams.minAreaUsingFlip = 100; %



algoParams=ReadYaml('/home/pdaude/Projet_Python/CREAM_PDFF/algoParams/B0NICE_algoParams.yml');

modelParams =ReadYaml('/home/pdaude/Projet_Python/CREAM_PDFF/modelParams/CustommodelParams.yml');
% Implemented in BONICE
%algoParams.species(2).frequency =[-242.7060 -217.1580 -166.0620 -123.9078 -24.9093 38.3220];
%algoParams.species(2).frequency  = algoParams.species(2).frequency *(3/1.5);
%img=load('Data_test/99.mat');

[species ,FWspectrum]= setupModelParams(modelParams);
algoParams.gyro=FWspectrum.gyro;
algoParams.species =species;

img=load('/home/pdaude/Projet_Python/CREAM_PDFF/fw_i3cm1i_3pluspoint_berglund_QPBO/test_cases/04.mat');
img.imDataParams.images = single(img.imDataParams.images);% Could not be integer

[params,sse,YM]=B0NICE_main(img.imDataParams,algoParams);
save('output_B0NICE.mat','params','sse','YM'); 
