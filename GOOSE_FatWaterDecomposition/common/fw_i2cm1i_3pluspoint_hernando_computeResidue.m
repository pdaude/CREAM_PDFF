% Objective: to calculate the residuel using the same formula as Hernando,
% ...........in order to compare the result from graph search algr.


%% Function name: fw_i2cm1i_3pluspoint_hernando
%%
%% Description: Fat-water separation using regularized fieldmap formulation and graph cut solution. 
%%
%% Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large 
%% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.
%% 
%% Some properties:
%%   - Image-space
%%   - 2 species (water-fat)
%%   - Complex-fitting
%%   - Multi-peak fat (pre-calibrated)
%%   - Single-R2*
%%   - Independent water/fat phase
%%   - Requires 3+ echoes at arbitrary echo times (some choices are much better than others! see NSA...)
%%
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%%   - imDataParams.TE: echo times (in seconds)
%%   - imDataParams.FieldStrength: (in Tesla)
%%
%%   - algoParams.species(ii).name = name of species ii (string)
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%%   Example
%%      - algoParams.species(1).name = 'water' % Water
%%      - algoParams.species(1).frequency = [0] 
%%      - algoParams.species(1).relAmps = [1]   
%%      - algoParams.species(2).name = 'fat' % Fat
%%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%% 
%%   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
%%   - algoParams.range_r2star = [0 0]; % Range of R2* values
%%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
%%   - algoParams.range_fm = [-400 400]; % Range of field map values
%%   - algoParams.NUM_FMS = 301; % Number of field map values to discretize
%%   - algoParams.NUM_ITERS = 40; % Number of graph cut iterations
%%   - algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
%%   - algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
%%   - algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
%%   - algoParams.lambda = 0.05; % Regularization parameter
%%   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
%%   - algoParams.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)  
%%   - algoParams.residual: in case we pre-computed the fit residual (mostly for testing) 
%%
%% Output: structure outParams
%%   - outParams.species(ii).name: name of the species (taken from algoParams)
%%   - outParams.species(ii).amps: estimated water/fat images, size [nx,ny,ncoils] 
%%   - outParams.r2starmap: R2* map (in s^{-1}, size [nx,ny])
%%   - outParams.fieldmap: field map (in Hz, size [nx,ny])
%%
%%
%% Author: Diego Hernando
%% Date created: August 5, 2011
%% Date last modified: November 10, 2011

function [residual, r2starmap] = fw_i2cm1i_3pluspoint_hernando_computeResidue( imDataParams, algoParams )


DEBUG = 0;

% Check validity of params, and set default algorithm parameters if not provided
[validParams,algoParams] = checkParamsAndSetDefaults( imDataParams,algoParams );
if validParams==0
  disp(['Exiting -- data not processed']);
  outParams = [];
  return;
end

% If precession is clockwise (positive fat frequency) simply conjugate data
if imDataParams.PrecessionIsClockwise <= 0 
  imDataParams.images = conj(imDataParams.images);
  imDataParams.PrecessionIsClockwise = 1;
end

% Check spatial subsampling option (speedup ~ quadratic SUBSAMPLE parameter)
SUBSAMPLE = algoParams.SUBSAMPLE;
if SUBSAMPLE > 1
  images0 = imDataParams.images;
  START = round(SUBSAMPLE/2);
  [sx,sy] = size(images0(:,:,1,1,1));
  allX = 1:sx;
  allY = 1:sy;
  subX = START:SUBSAMPLE:sx;
  subY = START:SUBSAMPLE:sy;
  imDataParams.images = images0(subX,subY,:,:,:);
end

% Regularization parameter
lambda=algoParams.lambda;

% Spatially-varying regularization.  The LMAP_POWER applies to the
% sqrt of the curvature of the residual, and LMAP_POWER=2 yields
% approximately uniform resolution.
LMAP_POWER = algoParams.LMAP_POWER;
  
  
% LMAP_EXTRA: Extra flexibility for including prior knowledge into
% regularization. For instance, it can be used to add more smoothing
% to noise regions (by adding, eg a constant LMAP_EXTRA), or even to
% add spatially-varying smoothing as a function of distance to
% isocenter...
LMAP_EXTRA = algoParams.LMAP_EXTRA;

% Finish off with some optimization transfer -- to remove discretization
DO_OT = algoParams.DO_OT;

% Let's get the residual. If it's not already in the params, compute it
try 
  % Grab the residual from the params structure
  residual = algoParams.residual;
catch
  % Check for uniform TE spacings
  dTE = diff(imDataParams.TE);
  
  try 
    TRY_PERIODIC_RESIDUAL = algoParams.TRY_PERIODIC_RESIDUAL;
  catch
    TRY_PERIODIC_RESIDUAL=1;
  end
  
  if TRY_PERIODIC_RESIDUAL==1 && sum(abs(dTE - dTE(1)))<1e-6 % If we have uniform TE spacing
    UNIFORM_TEs = 1;
  else
    UNIFORM_TEs = 0;
  end    


  if DEBUG == 1
    UNIFORM_TEs
  end

  % Compute the residual
  if UNIFORM_TEs == 1 % TEST DH* 090801
    % Find out the period in the residual (Assuming uniformly spaced samples)
    dt = imDataParams.TE(2)-imDataParams.TE(1);
    period = abs(1/dt);
    NUM_FMS_ORIG = algoParams.NUM_FMS;
    range = diff(algoParams.range_fm);
    params.NUM_FMS = ceil(algoParams.NUM_FMS/range*period);
    params.range_fm = [0 period*(1-1/(algoParams.NUM_FMS))];
    [residual,r2starmap] = GOOSE_computeResidual( imDataParams, algoParams );
    num_periods = ceil(range/period/2);
    algoParams.NUM_FMS = 2*num_periods*algoParams.NUM_FMS;
    residual = repmat(residual,[2*num_periods 1 1]);
    algoParams.range_fm = [-num_periods*period (num_periods*period-period/NUM_FMS_ORIG)];
    disp("range_fm")
    disp(algoParams.range_fm)
    disp("Num_FMS")
    disp(algoParams.NUM_FMS)
  else
    % If not uniformly spaced TEs, get the residual for the whole range
      [residual, r2starmap] = GOOSE_computeResidual( imDataParams, algoParams );
  end

end


  
end




