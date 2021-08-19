function [params,sse] = hernando_main(imDataParams,algoParams)
%% Recon -- graph cut 
%% (Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large 
%% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.)
[nx,ny,nz,ncoils,ne]=size(imDataParams.images);
imDataParamsZ=imDataParams;
for z =1:nz
imDataParamsZ.images = double(imDataParams.images(:,:,z,:,:));
outParams = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParamsZ, algoParams );


params.B0(:,:,z)= outParams.fieldmap ; % B0 (Hz)
params.R2(:,:,z) = outParams.r2starmap;
params.F(:,:,z) = outParams.species(2).amps ;
params.W(:,:,z) = outParams.species(1).amps ;
end

params.FF = squeeze(100*abs(params.F)./(abs(params.F)+abs(params.W))); % FF (%)
[sse,YM] = calculate_residual(params,algoParams,imDataParams);

end