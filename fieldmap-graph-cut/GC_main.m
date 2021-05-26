function [params, sse] = GC_main(imDataParams,algoParams)
%% check validaty of parameters
algoParams = checkParamsAndSetDefaults_GANDALF(imDataParams, algoParams);

%% calculate VARPRO residual and extract minima 
VARPROparams = calculateMinimaDirect(imDataParams, algoParams);

%% perform core graph cut function that outputs field-map, water, fat and r2*-map
waterFatParams = GANDALF(imDataParams, algoParams, VARPROparams);

params.B0 =waterFatParams.fieldmap ; % B0 (Hz)
params.R2 = waterFatParams.r2starmap; % R2* (1/s)

%params.PH =  ; % initial phase (rad)
params.F = waterFatParams.fat;
params.W = waterFatParams.water;
params.FF = squeeze(100*abs(params.F)./(abs(params.F)+abs(params.W))); % FF (%)

[sse,YM]=calculate_residual(params,algoParams,imDataParams);
end
