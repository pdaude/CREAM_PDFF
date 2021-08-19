function [params,sse,outParams]= Berglund_main(imDataParams,algoParams)
algoParams.fat_R2s = zeros(size(algoParams.species(2).relAmps));
outParams = fw_i3cm1i_3pluspoint_berglund_QPBO(imDataParams,algoParams);
params.B0= outParams.fieldmap ; % B0 (Hz)
params.R2 = outParams.r2starmap; %(s-1)
params.W = outParams.species(1).amps;
params.F = outParams.species(2).amps;
abs_water = abs(params.W);
abs_fat = abs(params.F);
FF = 100*(abs_fat)./(abs_fat + abs_water);
%FF(isnan(FF)) = 0;
params.FF=FF;
[sse,YM] = calculate_residual(params,algoParams,imDataParams);
end
