function [outParams] = GraphSearch(dataname, z)

%% load data
load(dataname);
imDataParams.images = double(imDataParams.images(:,:,z,:,:));
x = size(imDataParams.images, 1);
y = size(imDataParams.images, 2);
T = imDataParams.TE;
Numlayers = 100; % number of layers in the graph
algoParams.gyro=42.57747892;
gyro = algoParams.gyro%42.58;
freqRange = gyro*imDataParams.FieldStrength*8; %search range is [-8ppm, 8ppm];
t = linspace(-freqRange, freqRange, Numlayers);
gridspacing = t(2)-t(1);

%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
algoParams.gyro=42.57747892;

% Algorithm-specific parameters
algoParams.gridsize = round(gridspacing);
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 30]; % Range of R2* values
algoParams.NUM_R2STARS = 16; % Numbre of R2* values for quantization
algoParams.range_fm = [-freqRange freqRange]; % Range of field map values
algoParams.NUM_FMS = Numlayers; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 1; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;



%% compute VARPRO residue using hernando's code

[residual, r2starMap] = fw_i2cm1i_3pluspoint_hernando_computeResidue( imDataParams, algoParams );
residue = permute(residual,[2 3 1]);

%% graph search for fieldmap
% Header file to plug into surface detection

Header.Dataname = dataname;
Header.Delta = (T(end)-T(1))/(length(T)-1);
Header.freqRange = freqRange;
Header.fatpeaksNo = 6;

Header.stepsize = gridspacing;
Header.NumLayer = Numlayers;
Header.Z = z;
Header.Date = 203;
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
        fieldmap_rough(x,y) = -freqRange + optimal_surface(x,y)*gridspacing;
    end
end
[ fieldmapFine ] = findRangeMin( imDataParams, algoParams, fieldmap_rough, optimal_r2starmap );
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
  imDataParams.images = images0;
end
algoParams.lmap = lmap;

% Now take the field map fm and get the rest of the estimates
for ka=1:size(imDataParams.images,6) 
  
  curParams = imDataParams;
  curParams.images = imDataParams.images(:,:,:,:,:,ka);

  
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
  imDataParams.fmGC = fm;

  algoParams.OT_ITERS = 10;
  algoParams.fieldmap = fm;
  algoParams.r2starmap = r2starmap;
  algoParams.lambdamap = sqrt(algoParams.lambda*lmap);
  
  outParams = GOOSE_fw_i2cm0i_3plusploint_hernando_optimtransfer( imDataParams, algoParams  );      
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
end
