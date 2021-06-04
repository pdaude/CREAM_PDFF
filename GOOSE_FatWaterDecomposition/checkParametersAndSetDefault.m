function [algoParams] = checkParametersAndSetDefault(algoParams,Tesla)
%--------------
%set up default algoParams
% Algorithm-specific parameters


defaultalgoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
defaultalgoParams.range_r2star = [0 30]; % Range of R2* values
defaultalgoParams.NUM_R2STARS = 16; % Numbre of R2* values for quantization
defaultalgoParams.range_fm_ppm = [-8 8]; % Range of field map values in ppm (min / max)
defaultalgoParams.NUM_FMS = 100; % Number of field map values to discretize
defaultalgoParams.NUM_ITERS = 40; % Number of graph cut iterations
defaultalgoParams.SUBSAMPLE = 1; % Spatial subsampling for field map estimation (for speed)
defaultalgoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
defaultalgoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
defaultalgoParams.lambda = 0.05; % Regularization parameter
defaultalgoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
defaultalgoParams.TRY_PERIODIC_RESIDUAL = 0;
defaultalgoParams.gyro=42.57747892;
 
defaultfields=fieldnames(defaultalgoParams);
 
 for n= 1:numel(defaultfields)
     strfield=string(defaultfields(n));
     if ~isfield(algoParams,defaultfields(n))
         algoParams.(strfield)=defaultalgoParams.(strfield);
     end
 end
freqmin= algoParams.gyro*Tesla*algoParams.range_fm_ppm(1) ; %in Hz
freqmax= algoParams.gyro*Tesla*algoParams.range_fm_ppm(2) ; %in Hz
t = linspace(freqmin, freqmax, algoParams.NUM_FMS);
gridspacing = t(2)-t(1);
algoParams.range_fm = [freqmin freqmax];
algoParams.gridspacing = gridspacing;
algoParams.gridsize = round(gridspacing);
end