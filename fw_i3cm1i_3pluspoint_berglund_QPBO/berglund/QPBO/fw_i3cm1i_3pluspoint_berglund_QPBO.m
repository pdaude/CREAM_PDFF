%% Function name: fw_i3cm1i_3pluspoint_berglund_QPBO
%%
%% Description: Fat-water separation from three plus complex echoes with uniform
%%              echo time spacing, using a whole image optimization algorithm described in:
%%
%% Berglund J, Kullberg J. Three-dimensional water/fat separation and T2* estimation based on whole-image
%% optimization--application in breathhold liver imaging at 1.5 T. Magn Reson Med. 2012, Jun;67(6):1684-93.
%% 
%% Some properties:
%%   - Image-space
%%   - 2 species (water-fat)
%%   - Complex-fitting
%%   - Multi-peak fat (pre-calibrated)
%%   - R2* (given >3 echoes)
%%   - Independent water/fat phase
%%   - Requires 3 or more uniformly spaced echoes
%%
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,nz,ncoils,nTE]
%%   - imDataParams.TE: echo times (in seconds)
%%   - imDataParams.FieldStrength: (in Tesla)
%%   - imDataParams.voxelSize: (mm x mm x mm)
%%   - imDataParams.mask: image mask, array of size[nx,ny,nz]
%%   - imDataParams.PrecessionIsClockwise: ==1 will not conjugate, ~=1 will conjugate
%%
%%   - algoParams.species(ii).name = name of species ii (string)
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude of each peak within species ii
%%
%%   Example
%%   - algoParams.species(1).name = 'water';
%%   - algoParams.species(1).frequency = 4.70;
%%   - algoParams.species(1).relAmps = 1;
%%   - algoParams.species(2).name = 'fat';
%%   - algoParams.species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29]; % 9-peak model
%%   - algoParams.species(2).relAmps   = [  88,  642,   58,   62,   58,    6,   39,   10,   37]; % Hamilton G, et al. NMR Biomed. 24(7):784-90, 2011. PMID: 21834002
%% 
%%   - algoParams.decoupled_estimation = true; % flag for decoupled R2 estimation
%%   - algoParams.Fibonacci_search = true; % flag for Fibonacci search
%%   - algoParams.B0_smooth_in_stack_direction = false; % flag for B0 smooth in stack direction
%%   - algoParams.multigrid = true; % flag for multi-level resolution pyramid
%%   - algoParams.estimate_R2 = true; % flag to estimate R2star
%%   - algoParams.verbose = false; % flag for verbose status messages (default false)
%%   - algoParams.process_in_3D = true; % flag to process in 3D (default true)
%%   - algoParams.R2star_calibration = false; % flag to perform R2* calibration (default false)
%%   - algoParams.ICM_iterations = 2; % ICM iterations
%%   - algoParams.num_B0_labels = 100; % number of discretized B0 values
%%   - algoParams.mu = 10; % regularization parameter
%%   - algoParams.R2_stepsize = 1; % R2 stepsize in s^-1
%%   - algoParams.max_R2 = 120; % maximum R2 in s^-1
%%   - algoParams.max_label_change = 0.1; % 
%%   - algoParams.fine_R2_stepsize = 1.0; % Fine stepsize for discretization of R2(*) [s^-1] (used in decoupled fine-tuning step)
%%   - algoParams.coarse_R2_stepsize = 10.0; % Coarse stepsize for discretization of R2(*) [s^-1] (used for joint estimation step, value 0.0 --> Decoupled estimation only
%%   - algoParams.water_R2 = 0.0; % Water R2 [sec-1]
%%   - algoParams.fat_R2s = zeros(1,9); % fat peak R2s [s^-1]
%%   - algoParams.R2star_calibration_max = 800; % max R2* to use for calibration [s^-1] (default 800)
%%   - algoParams.R2star_calibration_cdf_threshold; % threshold for R2* calibration cumulative density function [0,1] (default 0.98)
%%
%% Output: structure outParams
%%   - outParams.species(ii).name: name of the species (taken from algoParams)
%%   - outParams.species(ii).amps: estimated water/fat images, size [nx,ny,nz] 
%%   - outParams.fieldmap: field map (in Hz, size [nx,ny,nz]), fieldmap is NOT unwrapped
%%   - outParams.r2starmap: R2* map (in s^{-1}, size [nx,ny,nz])
%%   - outParams.status: QPBO status map (size [nx,ny,nz])
%%   - outParams.residual: residual error map (size [nx,ny,nz])
%%   - outParams.mask: mask (size [nx,ny,nz])
%%   - outParams.execution_time_seconds
%%
%% Author: E. Brian Welch
%% Date created: January 23, 2013
%% Date last modified: May 23, 2014

function outParams = fw_i3cm1i_3pluspoint_berglund_QPBO( imDataParams, algoParams ),
time_start = tic;
outParams = [];

%% check validity of params, and set default algorithm parameters if not provided
[validParams, imDataParams, algoParams] = checkParamsAndSetDefaults_QPBO( imDataParams, algoParams );
if validParams==0
  disp(['Exiting -- data not processed']);
  return;
end

%% Build a signal input structure for DixonQPBO
DixonQPBO_input.TE = imDataParams.TE(imDataParams.idx_TEeven);
DixonQPBO_input.field_strength = imDataParams.FieldStrength;
DixonQPBO_input.voxelSize = imDataParams.voxelSize;

%% put water in species(1) and fat in species(2)
if strcmpi(algoParams.species(2).name,'water'),
    DixonQPBO_input.water_chemical_shift  = algoParams.species(2).frequency(:);
    DixonQPBO_input.fat_chemical_shifts   = algoParams.species(1).frequency(:);
    DixonQPBO_input.fat_amps              = algoParams.species(1).relAmps(:);
else
    DixonQPBO_input.water_chemical_shift  = algoParams.species(1).frequency(:);
    DixonQPBO_input.fat_chemical_shifts   = algoParams.species(2).frequency(:);
    DixonQPBO_input.fat_amps              = algoParams.species(2).relAmps(:);    
end
DixonQPBO_input.fat_R2s = algoParams.fat_R2s;

DixonQPBO_input.verbose = algoParams.verbose;
DixonQPBO_input.decoupled_estimation = algoParams.decoupled_estimation;
DixonQPBO_input.Fibonacci_search = algoParams.Fibonacci_search;
DixonQPBO_input.B0_smooth_in_stack_direction = algoParams.B0_smooth_in_stack_direction;
DixonQPBO_input.multigrid = algoParams.multigrid;

if ( (algoParams.estimate_R2) && (imDataParams.nTEeven>=4) ),
    DixonQPBO_input.estimate_R2 = true;
else
    DixonQPBO_input.estimate_R2 = false;
end

DixonQPBO_input.num_B0_labels = algoParams.num_B0_labels;
DixonQPBO_input.ICM_iterations = algoParams.ICM_iterations;
DixonQPBO_input.mu = algoParams.mu;
DixonQPBO_input.R2_stepsize = algoParams.R2_stepsize;
DixonQPBO_input.max_R2 = algoParams.max_R2;
DixonQPBO_input.max_label_change = algoParams.max_label_change;
DixonQPBO_input.fine_R2_stepsize = algoParams.fine_R2_stepsize;
DixonQPBO_input.coarse_R2_stepsize = algoParams.coarse_R2_stepsize;
DixonQPBO_input.water_R2 = algoParams.water_R2;
DixonQPBO_input.InplaneOverThroughplaneVoxelsize = imDataParams.voxelSize(1) / imDataParams.voxelSize(3);
DixonQPBO_input.te_used = imDataParams.TE(imDataParams.idx_TEeven(1));
DixonQPBO_input.dte_used = imDataParams.dTEeven;
DixonQPBO_input.use_num_echos = imDataParams.nTEeven;

%% dimensions
[nx ny nz ncoils nTE] = size(imDataParams.images);

%% detect crop indices
crop_cushion = 1;
idx_crop_x = {};
idx_crop_y = {};
for sl=1:nz,
    [idx_x idx_y] = find( imDataParams.mask(:,:,sl) );
    if length(idx_x)>0,
        idx_crop_x_start = max( 1,min(idx_x(:))-crop_cushion);
        idx_crop_x_stop  = min(nx,max(idx_x(:))+crop_cushion);
        idx_crop_y_start = max( 1,min(idx_y(:))-crop_cushion);
        idx_crop_y_stop  = min(ny,max(idx_y(:))+crop_cushion);
        idx_crop_x{sl} = [idx_crop_x_start:idx_crop_x_stop];
        idx_crop_y{sl} = [idx_crop_y_start:idx_crop_y_stop];
    else
        idx_crop_x{sl} = [];
        idx_crop_y{sl} = [];
    end
end

%% Call DixonQPBO via DixonApp
if exist('DixonApp.m')~=2,
    mfilename('fullpath')
    [pathstr,name,ext] = fileparts(mfilename('fullpath'));
    addpath( sprintf('%s/DixonApp/', pathstr) );
end

outParams.fieldmap = zeros([nx ny nz]);
outParams.r2starmap = zeros([nx ny nz]);
outParams.status = zeros([nx ny nz]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM R2* CALIBRATION IF ACTIVATED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( (algoParams.R2star_calibration) && (DixonQPBO_input.estimate_R2) ),
    
    if algoParams.verbose,
        disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: Performing R2* calibration with initial max R2* = %.2f', algoParams.R2star_calibration_max ) );
    end
    DixonQPBO_input.max_R2 = algoParams.R2star_calibration_max;    
    r2starmap_calibration = zeros([nx ny nz]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Call DixonQPBO via DixonApp
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if algoParams.process_in_3D,

        crop_x_3D_start = Inf;
        crop_x_3D_stop = 0;
        crop_y_3D_start = Inf;
        crop_y_3D_stop = 0;

        for sl=1:nz,
            if ( min(idx_crop_x{sl})<crop_x_3D_start ),
                crop_x_3D_start = min(idx_crop_x{sl});
            end
            if ( max(idx_crop_x{sl})>crop_x_3D_stop ),
                crop_x_3D_stop = max(idx_crop_x{sl});
            end
            if ( min(idx_crop_y{sl})<crop_y_3D_start ),
                crop_y_3D_start = min(idx_crop_y{sl});
            end
            if ( max(idx_crop_y{sl})>crop_y_3D_stop ),
                crop_y_3D_stop = max(idx_crop_y{sl});
            end            
        end

        idx_crop_x_3D = [crop_x_3D_start:crop_x_3D_stop];
        idx_crop_y_3D = [crop_y_3D_start:crop_y_3D_stop];

        if algoParams.verbose,
            disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: R2* Calibration: Cropping 3D to [ %d:%d , %d:%d , %d] (%.2f%%)', crop_x_3D_start, crop_x_3D_stop, crop_y_3D_start, crop_y_3D_stop, nz, 100*(length(idx_crop_x_3D)*length(idx_crop_y_3D)*nz)/(nx*ny*nz) ) );
        end

        DixonQPBO_input.images = imDataParams.images(idx_crop_x_3D,idx_crop_y_3D,:,1,imDataParams.idx_TEeven);
        DixonQPBO_output = DixonApp(DixonQPBO_input);

        r2starmap_calibration(idx_crop_x_3D,idx_crop_y_3D,:) = DixonQPBO_output.r2starmap;
        
    else

        for sl=1:nz,

            if algoParams.verbose,
                disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: R2* Calibration: Cropping slice %d to save %.2f%%', sl, 100*(1-length(idx_crop_x{sl})*length(idx_crop_y{sl})/(nx*ny)) ) );
            end

            DixonQPBO_input.images = imDataParams.images(idx_crop_x{sl},idx_crop_y{sl},sl,1,imDataParams.idx_TEeven);
            DixonQPBO_output_slice = DixonApp(DixonQPBO_input);

            r2starmap_calibration(idx_crop_x{sl},idx_crop_y{sl},sl) = DixonQPBO_output_slice.r2starmap;      
        end

    end
    
    % examine R2* histogram
    idx_mask_R2cal = find(imDataParams.mask(:)==1);
    numel_idx_mask_R2cal = numel(idx_mask_R2cal);
    num_bins = 1000;
    bins = linspace(0,algoParams.R2star_calibration_max,num_bins);
	[bin_counts , bin_centers] = hist(r2starmap_calibration(idx_mask_R2cal), bins );
    r2starmap_calibration_bin_cumsum_percent = (cumsum(bin_counts)/numel_idx_mask_R2cal);
    outParams.bin_cumsum_percent = r2starmap_calibration_bin_cumsum_percent;
	outParams.R2star_calibration_cdf_threshold = algoParams.R2star_calibration_cdf_threshold;
	idx_cumsum_above_threshold = find( r2starmap_calibration_bin_cumsum_percent >= outParams.R2star_calibration_cdf_threshold);
            
    % check to see if cdf_threshold should be relaxed
    while ( (idx_cumsum_above_threshold(1) == num_bins) && (outParams.R2star_calibration_cdf_threshold>0) ),
        outParams.R2star_calibration_cdf_threshold = outParams.R2star_calibration_cdf_threshold - 0.0150;
        idx_cumsum_above_threshold = find( r2starmap_calibration_bin_cumsum_percent >= outParams.R2star_calibration_cdf_threshold);
    end
    outParams.GLOBAL_maxR2  = ceil( bin_centers(idx_cumsum_above_threshold(1)) );
            
    % default to something close to reasonable if this cdf method completely bombs
    if (idx_cumsum_above_threshold(1) == num_bins),
        if (imDataParams.FieldStrength>1.6),
            outParams.GLOBAL_maxR2 = 350;
        else
            outParams.GLOBAL_maxR2 = 200;
        end
        warning( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: R2* Calibration: could not find a GLOBAL_maxR2 less than calibration_maxR2 ... using GLOBAL_maxR2  = %.2f', outParams.GLOBAL_maxR2));
    else
        if algoParams.verbose,
            disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: R2* Calibration: successfully satisified R2star_calibration_cdf_threshold = %.2f', outParams.R2star_calibration_cdf_threshold) );
        end
    end
    
    if algoParams.verbose,
        disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: R2* Calibration: GLOBAL_maxR2  = %.2f', outParams.GLOBAL_maxR2) );
    end
    
    % Update DixonQPBO_input with calibrated max R2* value
    DixonQPBO_input.max_R2 = outParams.GLOBAL_maxR2;
    
else
    if algoParams.verbose,
        disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: Skipping R2* calibration') );
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call DixonQPBO via DixonApp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if algoParams.process_in_3D,
    
    crop_x_3D_start = Inf;
    crop_x_3D_stop = 0;
    crop_y_3D_start = Inf;
    crop_y_3D_stop = 0;

    for sl=1:nz,
        if ( min(idx_crop_x{sl})<crop_x_3D_start ),
            crop_x_3D_start = min(idx_crop_x{sl});
        end
        if ( max(idx_crop_x{sl})>crop_x_3D_stop ),
            crop_x_3D_stop = max(idx_crop_x{sl});
        end
        if ( min(idx_crop_y{sl})<crop_y_3D_start ),
            crop_y_3D_start = min(idx_crop_y{sl});
        end
        if ( max(idx_crop_y{sl})>crop_y_3D_stop ),
            crop_y_3D_stop = max(idx_crop_y{sl});
        end            
    end

    idx_crop_x_3D = [crop_x_3D_start:crop_x_3D_stop];
    idx_crop_y_3D = [crop_y_3D_start:crop_y_3D_stop];

    if algoParams.verbose,
        disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: Cropping 3D to [ %d:%d , %d:%d , %d] (%.2f%%)', crop_x_3D_start, crop_x_3D_stop, crop_y_3D_start, crop_y_3D_stop, nz, 100*(length(idx_crop_x_3D)*length(idx_crop_y_3D)*nz)/(nx*ny*nz) ) );
    end
    
    DixonQPBO_input.images = imDataParams.images(idx_crop_x_3D,idx_crop_y_3D,:,1,imDataParams.idx_TEeven);
    DixonQPBO_output = DixonApp(DixonQPBO_input);
        
    outParams.fieldmap(idx_crop_x_3D,idx_crop_y_3D,:) = DixonQPBO_output.fieldmap * ( pi / 180.0 ) / (2*pi*imDataParams.dTEeven);
    outParams.r2starmap(idx_crop_x_3D,idx_crop_y_3D,:) = DixonQPBO_output.r2starmap;
    outParams.status(idx_crop_x_3D,idx_crop_y_3D,:) = DixonQPBO_output.status;
    
else
    
    for sl=1:nz,
        
        if algoParams.verbose,
            disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: Cropping slice %d to save %.2f%%', sl, 100*(1-length(idx_crop_x{sl})*length(idx_crop_y{sl})/(nx*ny)) ) );
        end
        
        DixonQPBO_input.images = imDataParams.images(idx_crop_x{sl},idx_crop_y{sl},sl,1,imDataParams.idx_TEeven);
        DixonQPBO_output_slice = DixonApp(DixonQPBO_input);
        
        outParams.fieldmap(idx_crop_x{sl},idx_crop_y{sl},sl) = DixonQPBO_output_slice.fieldmap * ( pi / 180.0 ) / (2*pi*imDataParams.dTEeven);
        outParams.r2starmap(idx_crop_x{sl},idx_crop_y{sl},sl) = DixonQPBO_output_slice.r2starmap;
        outParams.status(idx_crop_x{sl},idx_crop_y{sl},sl) = DixonQPBO_output_slice.status;        
    end
    
end

%% Decompose using fieldmap and r2starmap
algoParams_DECOMPOSE = algoParams;
try 
    ampW = algoParams_DECOMPOSE.species(1).relAmps;
catch
    ampW = 1.0;
end
gamma_Hz_per_Tesla = 42.577481e6;
deltaF = [0 ; gamma_Hz_per_Tesla/1e6*(algoParams_DECOMPOSE.species(2).frequency(:) - algoParams_DECOMPOSE.species(1).frequency(1))*(imDataParams.FieldStrength)];
relAmps = algoParams_DECOMPOSE.species(2).relAmps;
t = imDataParams.TE(:); % uses all echoes!
relAmps = reshape(relAmps,1,[]);  
        
B1 = zeros(nTE,2);
B = zeros(nTE,2);
outParams.species(1).amps = zeros([nx ny nz]);
outParams.species(2).amps = zeros([nx ny nz]);
outParams.residual = Inf*ones([nx ny nz]);

for n=1:nTE,
    B1(n,:) = [ampW*exp(1i*2*pi*deltaF(1)*t(n)),sum(relAmps(:).*exp(1i*2*pi*deltaF(2:end)*t(n)))];
end

for kz=1:nz,
    
    % cover entire crop boundary box
    idx_x = idx_crop_x{kz};
    idx_y = idx_crop_y{kz};

    for kx=idx_x,
        for ky=idx_y,

            s = reshape( imDataParams.images(kx,ky,kz,1,:), [nTE 1]);
            B(:,1) = B1(:,1).*exp(1i*2*pi*outParams.fieldmap(kx,ky,kz)*t(:) - outParams.r2starmap(kx,ky,kz)*t(:));
            B(:,2) = B1(:,2).*exp(1i*2*pi*outParams.fieldmap(kx,ky,kz)*t(:) - outParams.r2starmap(kx,ky,kz)*t(:));
            amps = B\s;
            outParams.species(1).amps(kx,ky,kz) = amps(1);
            outParams.species(2).amps(kx,ky,kz) = amps(2);
            outParams.residual(kx,ky,kz) = norm(s - B*amps,'fro');

        end % end ky
    end % end kx
end % end kz

%% Build outParams structure
outParams.species(1).name = algoParams.species(1).name;
outParams.species(2).name = algoParams.species(2).name;
if strcmpi(outParams.species(2).name,'water'),
    tmp = outParams.species(2).amps;
    outParams.species(2).amps = outParams.species(1).amps;
    outParams.species(1).amps = tmp;
end

outParams.mask = imDataParams.mask;
outParams.execution_time_seconds = toc(time_start);

end % end function fw_i3cm1i_3pluspoint_berglund_QPBO

function [validParams, imDataParams, algoParams] = checkParamsAndSetDefaults_QPBO( imDataParams, algoParams ),

validParams = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% imDataParams.verbose
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(algoParams,'verbose'),
        algoParams.verbose = false;
    else
        if ( (algoParams.verbose==1) || (algoParams.verbose==true) ),
            algoParams.verbose = true;
        else
            algoParams.verbose = false;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% imDataParams.images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(imDataParams,'images'),
        warning('fw_i3cm1i_3pluspoint_berglund_QPBO: input must contain field ''images'' with size [nx ny nz ncoils nTE]');
        return;
    else
        if ndims(imDataParams.images)~=5,
            warning('fw_i3cm1i_3pluspoint_berglund_QPBO: imDataParams.images should have dimensions [nx ny nz ncoils nTE]'); % also produces error for nTE=1
            return;
        else
            [nx ny nz ncoils nTE] = size(imDataParams.images);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% imDataParams.TE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(imDataParams,'TE'),
        warning('fw_i3cm1i_3pluspoint_berglund_QPBO: input must contain field ''TE'' with number of elements nTE');
        return;
    end
    if length(imDataParams.TE(:))~=nTE,
        warning('fw_i3cm1i_3pluspoint_berglund_QPBO; imDataParams.TE has length %d (expected %d)', length(imDataParams.TE(:)), nTE);
        return;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% imDataParams.FieldStrength
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(imDataParams,'FieldStrength'),
        imDataParams.FieldStrength = 1.5;
        warning('fw_i3cm1i_3pluspoint_berglund_QPBO: defaulting imDataParams.FieldStrength to 1.5');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% imDataParams.PrecessionIsClockwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(imDataParams,'PrecessionIsClockwise'),
        imDataParams.PrecessionIsClockwise = 1;
        warning('fw_i3cm1i_3pluspoint_berglund_QPBO: defaulting imDataParams.PrecessionIsClockwise to 1');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% imDataParams.mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(imDataParams,'mask'),
        imDataParams.mask = true([nx, ny, nz]);
        warning( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: defaulting imDataParams.mask to true([%d, %d, %d])', nx, ny, nz) );
    else
        imDataParams.mask = logical(imDataParams.mask);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% unhandled nTE values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nTE<3,
        warning( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: only %d echoes provided when at least 3 evenly spaced echoes are needed', nTE) );
        return;
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% combine coils
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (ncoils>1),
        imDataParams.images = coilCombine3D(imDataParams.images);
        ncoils = 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% account for precession direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if imDataParams.PrecessionIsClockwise ~= 1, 
        imDataParams.images = conj(imDataParams.images);
        imDataParams.PrecessionIsClockwise = 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% sort echo times (just to be sure)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imDataParams.TE = imDataParams.TE(:);
    [dummy,sorted_idx] = sort(imDataParams.TE);
    imDataParams.TE = imDataParams.TE(sorted_idx);
    imDataParams.images = imDataParams.images(:,:,:,:,sorted_idx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% average redundant echo times
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TE_usec_unique = unique( round(1e6*(imDataParams.TE(:))) );
    nTE_usec_unique = length(TE_usec_unique(:));
    if nTE_usec_unique~=nTE,
        if algoParams.verbose,
            disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: Found %d unique echo times out of %d provided echo times : Redundant echo times will be averaged', nTE_usec_unique, nTE) );
        end
        images_tmp = zeros([nx ny nz 1 nTE_usec_unique]);
        TE_tmp = zeros(nTE_usec_unique,1);
        for kTE=1:nTE_usec_unique,
            idx_TE = find(TE_usec_unique(kTE)==imDataParams.TE(:));
            images_tmp(:,:,:,1,kTE) = sum(imDataParams.images(:,:,:,1,idx_TE),5)/length(idx_TE(:));
            TE_tmp(kTE) = imDataParams.TE(idx_TE(1));
        end
        imDataParams.images = images_tmp;
        imDataParams.TE = TE_tmp;
        clear images_tmp TE_tmp;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% detect indices and number of evenly spaced echoes, idx_TEeven and nTEeven
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dTE_usec = round(1e6*diff(imDataParams.TE(:)));
    dTE_candidates = dTE_usec(:);
    for idx_dTE = 1:length(dTE_usec(:)),
        dTE_candidates = unique([dTE_candidates(:) ; cumsum(dTE_usec(idx_dTE:end)) ]);
    end
    dTE_candidates = unique(dTE_candidates(:));

    best_echoes = [];
    for idx_dTE_candidate = 1:length(dTE_candidates),
        target_dTE_usec = dTE_candidates(idx_dTE_candidate);
        ne = length(imDataParams.TE);

        % initialize dTE_groups
        dTE_groups = {};
        for idx=1:(ne-1),
            dTE_groups{idx}.dTE_usec = dTE_usec(idx);
            dTE_groups{idx}.echoes = [idx idx+1];
        end

        % combine dTE_groups
        for idx=1:(ne-1),
            if (dTE_groups{idx}.dTE_usec == target_dTE_usec),
                % do nothing
            elseif (dTE_groups{idx}.dTE_usec < target_dTE_usec)
                % combine with next group
                    if idx<(ne-1),
                        dTE_groups{idx+1}.dTE_usec = dTE_groups{idx}.dTE_usec + dTE_groups{idx+1}.dTE_usec;
                        dTE_groups{idx+1}.echoes = [dTE_groups{idx}.echoes(1:end-1) dTE_groups{idx+1}.echoes(2:end)];
                        dTE_groups{idx}.dTE_usec = NaN;
                        dTE_groups{idx}.echoes = [];
                    else
                        dTE_groups{idx}.dTE_usec = NaN;
                        dTE_groups{idx}.echoes = [];
                    end
            else
                % echo spacing is too large
                dTE_groups{idx}.dTE_usec = NaN;
                dTE_groups{idx}.echoes = [];
            end
        end

        % find usable TE set
        echoes_to_use = [];
        idx = 1;
        while idx <= (ne-1),
            if ( isnan(dTE_groups{idx}.dTE_usec) && (length(echoes_to_use)>0) ),
                break;
            end
            if ~isnan(dTE_groups{idx}.dTE_usec),
                echoes_to_use = [echoes_to_use dTE_groups{idx}.echoes];
            end
            idx = idx + 1;
        end

        echoes_to_use = sort( unique(echoes_to_use) );

        % test if this is better (longer) echo set
        if length(echoes_to_use) > length(best_echoes),
            best_echoes = echoes_to_use;
        end
    end

    imDataParams.idx_TEeven = best_echoes;
    imDataParams.nTEeven = length(imDataParams.idx_TEeven);

    % double check that detect evenly space echoes are truly evenly spaced
    dTEeven_usec = round(1e6*diff(imDataParams.TE(imDataParams.idx_TEeven)));
    if all(dTEeven_usec==dTEeven_usec(1)),
        imDataParams.dTEeven = imDataParams.TE(imDataParams.idx_TEeven(2)) - imDataParams.TE(imDataParams.idx_TEeven(1));
        if algoParams.verbose,
            disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: Available evenly space echoes confirmed with delta TE = %.6f seconds', imDataParams.dTEeven) );
        end
    end

    if imDataParams.nTEeven<3,
        warning('fw_i3cm1i_3pluspoint_berglund_QPBO: Could not find at least 3 evenly spaced echoes');
        return;
    else
        if algoParams.verbose,
            disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO: Found %d of %d evenly spaced echoes. idx_nTEeven = [ %s]', imDataParams.nTEeven, nTE, sprintf('%d ', imDataParams.idx_TEeven) ) );
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% check voxelsize
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if ~isfield(imDataParams,'voxelSize'),
        imDataParams.voxelSize = [1 1 1];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% QPBO default algorithm parameter values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(algoParams,'decoupled_estimation'),
        algoParams.decoupled_estimation = true;
    end
    if ~isfield(algoParams,'Fibonacci_search'),
        algoParams.Fibonacci_search = true;
    end
    if ~isfield(algoParams,'B0_smooth_in_stack_direction'),
        algoParams.B0_smooth_in_stack_direction = false;
    end
    if ~isfield(algoParams,'multigrid'),
        algoParams.multigrid = true;
    end
    if ~isfield(algoParams,'estimate_R2'),
        algoParams.estimate_R2 = true;
    end
    if ~isfield(algoParams,'verbose'),
        algoParams.verbose = false;
    end
    if ~isfield(algoParams,'process_in_3D'),
        algoParams.process_in_3D = true;
    end
    if ~isfield(algoParams,'R2star_calibration'),
        algoParams.R2star_calibration = false;
    end        
    if ~isfield(algoParams,'species'),
        algoParams.species(1).name = 'water';
        algoParams.species(1).frequency = 4.70;
        algoParams.species(1).relAmps = 1;
        algoParams.species(2).name = 'fat';
        algoParams.species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29]; % 9-peak model
        algoParams.species(2).relAmps   = [  88,  642,   58,   62,   58,    6,   39,   10,   37]; % Hamilton G, et al. NMR Biomed. 24(7):784-90, 2011. PMID: 21834002
    end
    algoParams.species(1).relAmps = algoParams.species(1).relAmps / sum( algoParams.species(1).relAmps(:) );
    algoParams.species(2).relAmps = algoParams.species(2).relAmps / sum( algoParams.species(2).relAmps(:) );
    
    if ~isfield(algoParams,'ICM_iterations'),
        algoParams.ICM_iterations = 2;
    end
    if ~isfield(algoParams,'num_B0_labels'),
        algoParams.num_B0_labels = 100;
    end
    if ~isfield(algoParams,'mu'),
        algoParams.mu = 10;
    end
    if ~isfield(algoParams,'R2_stepsize'),
        algoParams.R2_stepsize = 1;
    end
    if ~isfield(algoParams,'max_R2'),
        algoParams.max_R2 = 120;
    end
    if ~isfield(algoParams,'max_label_change'),
        algoParams.max_label_change = 0.1;
    end
    if ~isfield(algoParams,'fine_R2_stepsize'),
        algoParams.fine_R2_stepsize = 1.0;
    end
    if ~isfield(algoParams,'coarse_R2_stepsize'),
        algoParams.coarse_R2_stepsize = 10.0; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% QPBO check flag values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if algoParams.decoupled_estimation == 1,
        algoParams.decoupled_estimation = true;
    else
        algoParams.decoupled_estimation = false;
    end

    if algoParams.Fibonacci_search == 1,
        algoParams.Fibonacci_search = true;
    else
        algoParams.Fibonacci_search = false;
    end    

    if algoParams.B0_smooth_in_stack_direction == 1,
        algoParams.B0_smooth_in_stack_direction = true;
    else
        algoParams.B0_smooth_in_stack_direction = false;
    end    

    if algoParams.multigrid == 1,
        algoParams.multigrid = true;
    else
        algoParams.multigrid = false;
    end   

    if algoParams.estimate_R2 == 1,
        algoParams.estimate_R2 = true;
    else
        algoParams.estimate_R2 = false;
    end     
    
    if algoParams.verbose == 1,
        algoParams.verbose = true;
    else
        algoParams.verbose = false;
    end
    
	if algoParams.process_in_3D == 1,
        algoParams.process_in_3D = true;
    else
        algoParams.process_in_3D = false;
    end
    
    if algoParams.R2star_calibration == 1,
        algoParams.R2star_calibration = true;
    else
        algoParams.R2star_calibration = false;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% other default values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(algoParams,'water_R2'),
        algoParams.water_R2  = 0.0;
    end
    
    if ~isfield(algoParams,'fat_R2s'),
        if strcmpi(algoParams.species(2).name,'water'),
            algoParams.fat_R2s  = zeros(1,length(algoParams.species(1).frequency));
        else
            algoParams.fat_R2s  = zeros(1,length(algoParams.species(2).frequency));
        end
    end
    
    if ~isfield(algoParams,'R2star_calibration_max'),
        algoParams.R2star_calibration_max  = 800;
    end
    
    
    if ~isfield(algoParams,'R2star_calibration_cdf_threshold'),
        algoParams.R2star_calibration_cdf_threshold  = 0.98;
    end        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    validParams = 1;
    
end % end function checkParamsAndSetDefaults_QPBO

% Function: coilCombine
%
% Description: combine multi-coil image sequences
%
% Based on: Walsh DO, Gmitro AF, Marcellin MW. Adaptive reconstruction of
% phased array MR imagery. Magn Reson Med 2000;43:682-690
% 
% Parameters:
% im1: the multi-coil images (size [nx,ny,nz,ncoils,nTE])
%
% Returns: 
% im2: the coil-combined images (size [nx,ny,nz,1,nTE])
%
% Author: Diego Hernando
% Date created: August 13, 2011
% Date last modified: December 8, 2011

function im2 = coilCombine( im1 )

% Let's make the coil dimension the fourth one and the TE the third
im1 = permute(im1,[1 2 5 4 3]);

% Get image dimensions and set filter size
[sx,sy,N,C] = size(im1);
filtsize = 7;

% Initialize
im2 = zeros(sx,sy,1,1,N);
Rs = zeros(sx,sy,C,C);

% Get correlation matrices
for kc1=1:C
  for kc2=1:C
    for kn=1:N
      Rs(:,:,kc1,kc2) = Rs(:,:,kc1,kc2) + filter2(ones(filtsize),im1(:,:,kn,kc1).*conj(im1(:,:,kn,kc2)),'same');
    end
  end
end

% Compute and apply filter at each voxel
for kx=1:sx
  for ky=1:sy
% $$$     [U,S] = eig(squeeze(Rs(kx,ky,:,:)));
% $$$     s = diag(S);
% $$$     [maxval,maxind] = max(abs(s));
% $$$     myfilt = U(:,maxind);    
% $$$     im2(kx,ky,:) = myfilt'*reshape(squeeze(im1(kx,ky,:,:)).',[C N]);

    % Change suggested by Mark Bydder
    [U,S] = svd(squeeze(Rs(kx,ky,:,:)));
    myfilt = U(:,1); 
    im2(kx,ky,1,1,:) = myfilt'*reshape(squeeze(im1(kx,ky,:,:)).',[C N]);
  end
end

% In case the input data are single
if strcmp(class(im1),'single')
  im2 = single(im2);
end

end % end coilCombine

% Function: coilCombine3D
%
% Description: combine multi-coil image sequences
%
% Based on: Walsh DO, Gmitro AF, Marcellin MW. Adaptive reconstruction of
% phased array MR imagery. Magn Reson Med 2000;43:682-690
% 
% Parameters:
% im1: the multi-coil images (size [nx,ny,ncoils,nTE])
%
% Returns: 
% im2: the coil-combined images (size [nx,ny,nTE])
%
% Author: Diego Hernando
% Date created: August 13, 2011
% Date last modified: February 7, 2012

function im2 = coilCombine3D( im1 )


[sx,sy,sz,C,N] = size(im1);

% Maintain the data type (e.g., single, double) of the input data
ims = zeros([sx,sy,sz,1,N],class(im1));
for kz=1:sz
  im2(:,:,kz,1,:) = coilCombine(im1(:,:,kz,:,:));
end

end % end coilCombine3D