%% setup algoParams
clc
close all
clear all
algoParams=ReadYaml('/home/pdaude/Projet_Python/CREAM_PDFF/algoParams/Berglund_algoParams.yml');

%gamma_Hz_per_Tesla = 42.577481e6;
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 4.70;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29]; % 9-peak model
algoParams.species(2).relAmps   = [  88,  642,   58,   62,   58,    6,   39,   10,   37]; % Hamilton G, et al. NMR Biomed. 24(7):784-90, 2011. PMID: 21834002

img=load('/home/pdaude/Projet_Python/CREAM_PDFF/fw_i3cm1i_3pluspoint_berglund_QPBO/test_cases/01.mat');
[params,sse,outParams]=Berglund_main(img.imDataParams,algoParams);
save('output_Berglund.mat','params','sse','outParams'); 