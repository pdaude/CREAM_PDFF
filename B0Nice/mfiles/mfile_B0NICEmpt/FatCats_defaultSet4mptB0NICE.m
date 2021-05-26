function [algoParams] = FatCats_defaultSet4mptB0NICE(algoParams)
%--------------
%set up default algoParams
%for initial B0 mapping
defaultalgoParams.FILTsize = 0;
defaultalgoParams.th4unwrap = 0; %prctile slice by slice
defaultalgoParams.th4supp = 20;%prctile slice by slice
defaultalgoParams.th4STAS = 40;%prctile slice by slice
defaultalgoParams.unwrap_mode = 2; % 2 for 2D PUROR; 3 for 3D PUROR
% for magnitude fitting
defaultalgoParams.FF_trl = [0.0 0.15 0.25 0.5 0.75 0.85 1.0]; 
% defaultalgoParams.FF_trl = 0:0.1:1.0;
defaultalgoParams.R2star_pt = 0:1:300; % unit: Hz
% defaultalgoParams.FFth4magWater = 0.3;
%---------------------------------------------
% phase error correction
defaultalgoParams.flag_debug = 0;
defaultalgoParams.th4RegionDivision = pi/2;
defaultalgoParams.minAreaUsingPhaseGradient = 20; % 
defaultalgoParams.minAreaUsingFlip = 100; %
 
defaultfields=fieldnames(defaultalgoParams);
 
 for n= 1:numel(defaultfields)
     strfield=string(defaultfields(n));
     if ~isfield(algoParams,defaultfields(n))
         algoParams.(strfield)=defaultalgoParams.(strfield);
     end
 end
     
algoParams.CutOff = length(algoParams.TE_seq);

% % to remove the long TEs
% if  algoParams.B0_strength == 3
%      algoParams.CutOff = 5;
% else
%      algoParams.CutOff = length( defaultalgoParams.TE_seq);
end