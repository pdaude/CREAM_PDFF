function [algoParams] = FatCats_defaultSet4mptB0NICE(algoParams)
%--------------
%set up default algoParams
%for initial B0 mapping
    algoParams.FILTsize = 0;
    algoParams.th4unwrap = 0; %prctile slice by slice
    algoParams.th4supp = 20;%prctile slice by slice
    algoParams.th4STAS = 40;%prctile slice by slice
    algoParams.unwrap_mode = 2; % 2 for 2D PUROR; 3 for 3D PUROR
% for magnitude fitting
algoParams.FF_trl = [0.0 0.15 0.25 0.5 0.75 0.85 1.0]; 
%algoParams.FF_trl = 0:0.1:1.0;
algoParams.R2star_pt = 0:1:300; % unit: Hz
%algoParams.FFth4magWater = 0.3;
%---------------------------------------------
% phase error correction
    algoParams.flag_debug = 0;
    algoParams.th4RegionDivision = pi/2;
    algoParams.minAreaUsingPhaseGradient = 20; % 
    algoParams.minAreaUsingFlip = 100; %
% to remove the long TEs
if algoParams.B0_strength == 3
    algoParams.CutOff = 5;
else
    algoParams.CutOff = length(algoParams.TE_seq);
end