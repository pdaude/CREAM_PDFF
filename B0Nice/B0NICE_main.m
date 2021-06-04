function [params,sse,YM] = B0NICE_main(imDataParams,algoParams)

algoParams.complex_image = imDataParams.images;
algoParams.B0_strength = imDataParams.FieldStrength; % unit:T
algoParams.PrecessionIsClockwise = imDataParams.PrecessionIsClockwise;
algoParams.TE_seq=imDataParams.TE;% unit:s

%Transform fat frequency in ppm  to Hz and shiffted by water frequency
larmor=algoParams.B0_strength*algoParams.gyro;
algoParams.species(2).frequency = (algoParams.species(2).frequency-algoParams.species(1).frequency)*larmor;


[algoParams] = FatCats_defaultSet4mptB0NICE(algoParams);
[algoParams] = FatCats_VerifyInputData(algoParams);

%------------------------------
%step 1: initial B0 mapping
disp('Step 1 : initial B0 mapping');
algoParams.echo_selection = 0; % 0 auto, 1 manual
[algoParams] = FatCats_initialB0PhaseMapping_SelectEcho(algoParams);
%step 2: 
[algoParams] = FatCats_BuildComplexB0(algoParams);
%step 3: 
%generating masks used ofr phase unwrapping
[ algoParams] = initial_b0MAPPING(algoParams );
%step 4: generating mag_based fat/water mask
disp('Step 2 : generating mag_based fat/water mask');
%pixel-by-pixel magnitude fitting
% generate R2* maps and mag-based fat and water masks
[algoParams] = FatCats_MagFitting(algoParams);
%
%step 3: phase error correction 
disp('Step 3 : phase error correction');
[algoParams] = FatCats_PhaseErrorCorrection(algoParams);  
%#############################################
%Step 6: final B0 & fat-water separation
disp('Step 4 : final B0 & fat-water separation');
[outParams] = FatCats_BiasRemoval( algoParams );
params.B0 =outParams.fm ; % B0 (Hz)
params.R2 = outParams.R2star ; % R2* (1/s)
params.FF = outParams.FF*100; % FF (%)
%params.PH =  ; % initial phase (rad)
params.F = outParams.F;
params.W = outParams.W;
[sse,YM] = calculate_residual(params,algoParams,imDataParams);

end
