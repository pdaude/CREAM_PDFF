%% Test fw_i3cm1i_3pluspoint_berglund_QPBO
%clear all; close all; clc;

%% start clock
test_berglund_QPBO_start_time = tic;

%% detect run location details
[mfile_pathstr, mfile_name, mfile_ext] = fileparts(mfilename('fullpath'));

%% check for test_results folder
test_results_folder = sprintf('%s/test_results', mfile_pathstr);
if exist(test_results_folder)~=7,
    mkdir(test_results_folder);
end

%% start diary
diary_file_name = sprintf('%s/test_berglund_QPBO_diary.txt', test_results_folder);
warning off; delete(diary_file_name); warning on;
diary(diary_file_name);
diary on;

%% folder locations
case_folder = sprintf('%s/test_cases', mfile_pathstr);

%% setup algoParams
gamma_Hz_per_Tesla = 42.577481e6;
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 4.70;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29]; % 9-peak model
algoParams.species(2).relAmps   = [  88,  642,   58,   62,   58,    6,   39,   10,   37]; % Hamilton G, et al. NMR Biomed. 24(7):784-90, 2011. PMID: 21834002
algoParams.decoupled_estimation = true; % flag for decoupled R2 estimation
algoParams.Fibonacci_search = true; % flag for Fibonacci search
algoParams.B0_smooth_in_stack_direction = false; % flag for B0 smooth in stack direction
algoParams.multigrid = true; % flag for multi-level resolution pyramid
algoParams.estimate_R2 = true; % flag to estimate R2star
algoParams.verbose = true; % flag for verbose status messages (default false)
algoParams.process_in_3D = false; % flag to process in 3D (default true)
algoParams.R2star_calibration = true; % flag to perform R2* calibration (default false)
algoParams.ICM_iterations = 2; % ICM iterations
algoParams.num_B0_labels = 100; % number of discretized B0 values
algoParams.mu = 10; % regularization parameter
algoParams.R2_stepsize = 1; % R2 stepsize in s^-1
algoParams.max_R2 = 120; % maximum R2 in s^-1
algoParams.max_label_change = 0.1; % 
algoParams.fine_R2_stepsize = 1.0; % Fine stepsize for discretization of R2(*) [s^-1] (used in decoupled fine-tuning step)
algoParams.coarse_R2_stepsize = 10.0; % Coarse stepsize for discretization of R2(*) [s^-1] (used for joint estimation step, value 0.0 --> Decoupled estimation only
algoParams.water_R2 = 0.0; % Water R2 [sec-1]
algoParams.fat_R2s = zeros(1,9); % fat peak R2s [s^-1]
algoParams.R2star_calibration_max = 800; % max R2* to use for calibration [s^-1] (default 800)
algoParams.R2star_calibration_cdf_threshold = 0.98; %% threshold for R2* calibration cumulative density function [0,1] (default 0.98)

%% loop through cases
%for c=4, % Quick test
%for c=1:10, % Phase I cases
%for c=11:17, % Phase II cases
for c=4, % All cases
    
    clear imDataParams;
    case_matfilename = sprintf('%s/%02d.mat', case_folder, c);
    load(case_matfilename);
    [nx ny nz ncoils nTE] = size(imDataParams.images);
        
    outParams = fw_i3cm1i_3pluspoint_berglund_QPBO(imDataParams,algoParams);
    
    FW = cat(3, abs(outParams.species(2).amps), abs(outParams.species(1).amps) );
    FW = FW - min(FW(:));
    FW = FW / max(FW(:));
    F = FW(:,:,[1:nz]);
    W = FW(:,:,[1:nz]+nz);
    
    figure(1000+c);
    imagesc([F(:,:) ; W(:,:)]);
    axis image;
    colormap(gray);
    title(sprintf('CASE %02d',c));
    drawnow;
    
    FWpngfilename = sprintf('%s/FW_%02d_fw_i3cm1i_3pluspoint_berglund_QPBO.png', test_results_folder, c);
    imwrite([F(:,:) ; W(:,:)], FWpngfilename, 'PNG', 'BitDepth', 8, 'Author', 'fw_i3cm1i_3pluspoint_berglund_QPBO', 'Description', sprintf('2012 ISMRM Fat Water Challenge : Case %02d : %s', c, datestr(now)) );
    
    FSF = F./(F+W+eps);
    
    % to reduce noise bias, recalculate FSF in voxels where fat is not the dominant species
    idx_F_not_dominant = find( FSF(:)<0.5 );
    FSF(idx_F_not_dominant) = 1 - W(idx_F_not_dominant)./(F(idx_F_not_dominant)+W(idx_F_not_dominant)+eps);
    
    FSFpngfilename = sprintf('%s/FSF_%02d_fw_i3cm1i_3pluspoint_berglund_QPBO.png', test_results_folder, c);
    imwrite(256*FSF(:,:), rainbow, FSFpngfilename, 'PNG', 'BitDepth', 8, 'Author', 'fw_i3cm1i_3pluspoint_berglund_QPBO', 'Description', sprintf('2012 ISMRM Fat Water Challenge : Case %02d : %s', c, datestr(now)) );
    
    FSFpngfilename = sprintf('%s/FSF_%02d_fw_i3cm1i_3pluspoint_berglund_QPBO.png', test_results_folder, c);
    imwrite(256*FSF(:,:), rainbow, FSFpngfilename, 'PNG', 'BitDepth', 8, 'Author', 'fw_i3cm1i_3pluspoint_berglund_QPBO', 'Description', sprintf('2012 ISMRM Fat Water Challenge : Case %02d : %s', c, datestr(now)) );
    
    R2 = outParams.r2starmap;
    R2 = R2 / max(R2(:));
    R2pngfilename = sprintf('%s/R2_%02d_fw_i3cm1i_3pluspoint_berglund_QPBO.png', test_results_folder, c);
    imwrite(R2(:,:), R2pngfilename, 'PNG', 'BitDepth', 8, 'Author', 'fw_i3cm1i_3pluspoint_berglund_QPBO', 'Description', sprintf('2012 ISMRM Fat Water Challenge : Case %02d : %s', c, datestr(now)) );
    
    FM = outParams.fieldmap;
    FM = FM - min(FM(:));
    FM = FM / max(FM(:));
    FMpngfilename = sprintf('%s/FM_%02d_fw_i3cm1i_3pluspoint_berglund_QPBO.png', test_results_folder, c);
    imwrite(FM(:,:), FMpngfilename, 'PNG', 'BitDepth', 8, 'Author', 'fw_i3cm1i_3pluspoint_berglund_QPBO', 'Description', sprintf('2012 ISMRM Fat Water Challenge : Case %02d : %s', c, datestr(now)) );
    
    disp( sprintf('COMPLETED CASE %02d',c) );

end

%% completion time
disp(sprintf('Completed test_berglund_QPBO in %0.2f seconds', toc(test_berglund_QPBO_start_time) ));

%% stop diary
diary off;
