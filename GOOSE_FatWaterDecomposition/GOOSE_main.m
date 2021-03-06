function [params,sse] = GOOSE_main(imDataParams,algoParams)
[nx,ny,nz,ncoils,ne]=size(imDataParams.images);
imDataParamsZ=imDataParams;
[algoParams] = checkParametersAndSetDefault(algoParams,imDataParams.FieldStrength);
for z =1:nz
imDataParamsZ.images = double(imDataParams.images(:,:,z,:,:));

%% Set recon parameters
% General parameters


% Algorithm-specific parameters
% algoParams.gridsize = round(gridspacing);
% 
% algoParams.range_fm = [-freqRange freqRange]; % Range of field map values

%% compute VARPRO residue using hernando's code

[residual, r2starMap] = fw_i2cm1i_3pluspoint_hernando_computeResidue( imDataParamsZ, algoParams );
residue = permute(residual,[2 3 1]);

%% graph search for fieldmap
% Header file to plug into surface detection

%Header.Dataname = dataname;
% Header.Delta = (T(end)-T(1))/(length(T)-1);
% Header.freqRange = freqRange;
% Header.fatpeaksNo = 6;
% 
% Header.Z = z;
% Header.Date = 203;
Header.stepsize = algoParams.gridspacing;
Header.NumLayer = algoParams.NUM_FMS;
Header.Num2Star = algoParams.NUM_R2STARS;
Header.constraint =3;
Header.residueMin = min(residue(:));
Header.residueMax = max(residue(:));

%mex -I/nfs/s-iibi52/users/ccui1/Downloads/graphsearchwithKanglib/toolbox_for_submission detect_optimal_surface.cpp
% This is the core of the algorithm
[optimal_surface] = detect_optimal_surface(Header, residue);
optimal_r2starmap = zeros( size(optimal_surface,1), size(optimal_surface,2) );
fieldmap_rough = zeros( size(optimal_surface,1), size(optimal_surface,2) );
R = permute(squeeze(r2starMap),[2,3,1]);
for y = 1: size(optimal_surface,2)
    for x = 1:size(optimal_surface,1) 
        optimal_r2starmap(x,y) = R(x,y,optimal_surface(x,y));
        fieldmap_rough(x,y) = algoParams.range_fm(1) + optimal_surface(x,y)*algoParams.gridspacing;
    end
end
[ fieldmapFine ] = findRangeMin( imDataParamsZ, algoParams, fieldmap_rough, optimal_r2starmap );
% fieldmapFine is the fieldmap acquired from graph surface search

%% get w/f images
fms = linspace(algoParams.range_fm(1),algoParams.range_fm(2),algoParams.NUM_FMS);
dfm = fms(2)-fms(1);
lmap = GOOSE_getQuadraticApprox( residual, dfm );  
lmap = (sqrt(lmap)).^algoParams.LMAP_POWER;
lmap = lmap + mean(lmap(:))*algoParams.LMAP_EXTRA;

fm = fieldmapFine;
% If we have subsampled (for speed), let's interpolate the field map
if algoParams.SUBSAMPLE>1
  fmlowres = fm;
  [SUBX,SUBY] = meshgrid(subY(:),subX(:));
  [ALLX,ALLY] = meshgrid(allY(:),allX(:));
  fm = interp2(SUBX,SUBY,fmlowres,ALLX,ALLY,'*spline');
  lmap = interp2(SUBX,SUBY,lmap,ALLX,ALLY,'*spline');
  fm(isnan(fm)) = 0;
  imDataParamsZ.images = images0;
end
algoParams.lmap = lmap;

% Now take the field map fm and get the rest of the estimates
for ka=1:size(imDataParamsZ.images,6) 
  
  curParams = imDataParamsZ;
  curParams.images = imDataParamsZ.images(:,:,:,:,:,ka);

  
  if algoParams.range_r2star(2)>0
    % DH* 100422 use fine R2* discretization at this point 
    algoParams.NUM_R2STARS = round(algoParams.range_r2star(2)/2)+1; 
    r2starmap(:,:,ka) = GOOSE_estimateR2starGivenFieldmap( curParams, algoParams, fm );
  else
    r2starmap(:,:,ka) = zeros(size(fm));
  end
  
  if algoParams.DO_OT ~= 1
    % If no Optimization Transfer, just get the water/fat images
    amps = GOOSE_decomposeGivenFieldMapAndDampings( curParams,algoParams, fm,r2starmap(:,:,ka),r2starmap(:,:,ka) );
    waterimage = squeeze(amps(:,:,1,:));
    fatimage = squeeze(amps(:,:,2,:));
    w(:,:,:,ka) = waterimage;
    f(:,:,:,ka) = fatimage;
  end  
end

% If Optimization Transfer is requested, do it now
if algoParams.DO_OT == 1
  imDataParamsZ.fmGC = fm;

  algoParams.OT_ITERS = 10;
  algoParams.fieldmap = fm;
  algoParams.r2starmap = r2starmap;
  algoParams.lambdamap = sqrt(algoParams.lambda*lmap);
  
  outParams = GOOSE_fw_i2cm0i_3plusploint_hernando_optimtransfer( imDataParamsZ, algoParams  );      
  fm = outParams.fieldmap;
  

  % Now re-estimate the R2* map and water/fat images
  if algoParams.range_r2star(2)>0
    algoParams.NUM_R2STARS = round(algoParams.range_r2star(2)/2)+1; 
    r2starmap(:,:,ka) = GOOSE_estimateR2starGivenFieldmap( curParams,algoParams, fm );
  else
    r2starmap(:,:,ka) = zeros(size(fm));
  end
  % track r2starmap
  r2star(:,:,1) = optimal_r2starmap;
  r2star(:,:,2) = r2starmap(:,:,ka);
  
  % Now estimate the water/fat images
  amps = GOOSE_decomposeGivenFieldMapAndDampings( curParams,algoParams, fm,r2starmap(:,:,ka),r2starmap(:,:,ka) );
  waterimage = squeeze(amps(:,:,1,:));
  fatimage = squeeze(amps(:,:,2,:));
  outParams.water = waterimage;
  outParams.fat = fatimage;
  
 params.B0(:,:,z)= fm ; % B0 (Hz)
params.R2(:,:,z) = r2starmap(:,:,1);

params.F(:,:,z) = fatimage;
params.W(:,:,z) = waterimage;
end
end
params.FF = squeeze(100*abs(params.F)./(abs(params.F)+abs(params.W))); % FF (%)
[sse,YM] = calculate_residual(params,algoParams,imDataParams);
end
