clear all 
clc
%Create ground reference for wrong spectrum 
%Transform simuGt in wrong simuGT
gitPath='/home/pdaude/Projet_Python/CREAM_PDFF/';
imDataParamsFolder='simuImDataParams';
gtFolder='simuGT';
WgtFolder='WrongsimuGT';
modelParamsFolder='modelParams';
wrong_spectrum='HamiltonLiver2011';
simuPath=fullfile(gitPath,'simulation');
imDataParamsPath=fullfile(simuPath,imDataParamsFolder);
gtPath=fullfile(simuPath,gtFolder);
WgtPath=fullfile(simuPath,WgtFolder);

if not(isfolder(WgtPath))
    mkdir(WgtPath)
end

modelParamsPath=fullfile(gitPath,modelParamsFolder,join([char(wrong_spectrum),modelParamsFolder,'.yml']));
disp(modelParamsPath)
modelParams =ReadYaml(modelParamsPath);
[species ,FWspectrum]= setupModelParams(modelParams);
algoParams.gyro=FWspectrum.gyro;
algoParams.species=species;

list_gt = dir(fullfile(gtPath,'*.mat'));

for n=1:size(list_gt,1)
    path_gt=fullfile(list_gt(n).folder,list_gt(n).name);
    [~,gtfile,~]=fileparts(list_gt(n).name);
    splitted_gtfile=split(gtfile,'_');
    paramsExt=strjoin(splitted_gtfile(2:end),'_');
    simuImDataParamsPath=fullfile(imDataParamsPath,sprintf('%s_%s',imDataParamsFolder,paramsExt));
    load(path_gt,'params');
    load(simuImDataParamsPath,'imDataParams');
    [sse,~]=calculate_residual(params,algoParams,imDataParams);
    WsimuGTPath=fullfile(WgtPath,sprintf('Wrong%s',gtfile));
    save(WsimuGTPath,'params','sse','algoParams');
end 
