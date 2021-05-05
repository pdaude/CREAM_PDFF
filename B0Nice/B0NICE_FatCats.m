%
%% % % % % % % % % % % % % % % % % % % % % 
%       Junmin Liu (junminliumri@gmail.com)
%       University of Western Ontario, Canada
%       Dec 22, 2016
% % % % % % % % % % % % % % % % % % % % %
% modified by Junmin liu Dec 22, 2016
%  
%
%   Reference:
%   ---------
%   B0-NICE for B0 mapping
%   Junmin Liu and Maria Drangova,
%   "Method for B0 Off-Resonance Mapping by Non-Iterative Correction of
%   Phase-Errors (B0-NICE)
%   MRM. DOI 10.1002/mrm.25497
%
%   PUROR for phase unwrapping
%   Junmin Liu and Maria Drangova, 
%   "Intervention-based multidimensional phase unwrapping using recursive 
%   orthogonal referring",
%   Magnetic Resonance in Medicine, Volume 68(4):1303-1316, 2012
%   PUROR 3D and MEX files are also avalable upon request
%--------------
% Note that part of the code that related to refine the B0 and FF has been removed
% in order to make the frame work more clear but resulting in lower scores 
% for some of cases than that reported in MRM.25497.
%clc
clear all;close all;clc
addpath(genpath(pwd));
%
for index_case = 99:99 % example Junmin Liu et al., ISMRM 2014, P1615
%for index_case = 1:17 % will test all 17 cases from Challenge (www.ismrm.org/challenge/node/4)    
    case_str = num2str(index_case);
    if index_case < 10
    input_file = ['0',case_str,'.mat'];    
    output_FF_file = ['FF_','0',case_str,'_B0NICE.mat'];
    output_fm_file = ['fm_','0',case_str,'_B0NICE.mat'];
    else
    input_file = [case_str,'.mat'];    
    output_FF_file = ['FF_',case_str,'_B0NICE.mat'];
    output_fm_file = ['fm_',case_str,'_B0NICE.mat'];        
    end
    %
%-------------------------------------------------------------------
load(input_file);
if index_case == 99
imDataParams.images = single(imDataParams.images);
end
%
algoParams.complex_image = imDataParams.images;
algoParams.B0_strength = imDataParams.FieldStrength; % unit:T
algoParams.PrecessionIsClockwise = imDataParams.PrecessionIsClockwise;
algoParams.TE_seq = imDataParams.TE;% unit:s
%clear imDataParams;
%
[algoParams] = FatCats_defaultSet4mptB0NICE(algoParams);
[algoParams] = FatCats_VerifyInputData(algoParams);
algoParams.debug = 2; % 0 off all the figure; 1 on all; 2 only on the final
%------------------------------
%step 1: initial B0 mapping
disp('Step 1 : initial B0 mapping');
algoParams.echo_selection = 1; % 0 auto, 1 manual
algoParams.index_B0_suggested =[1 2 6;1 2 4;1 1 2;1 2 6;1 1 4;1 1 4;1 2 5;1 2 5;...
    1 1 6;1 2 5;1 1 3;1 1 4;1 2 4;1 1 3;1 1 4;1 2 5;2 1 3];
if index_case <= 17
disp('suggested echoes as : flag_B0 index_start index_end');
    index_B0 = algoParams.index_B0_suggested(index_case,:)
elseif index_case == 99
    disp('suggested echoes as : flag_B0 index_start index_end');
    index_B0 = [1 2 4]
end
%
[algoParams] = FatCats_initialB0PhaseMapping_SelectEcho(algoParams);
[algoParams] = FatCats_BuildComplexB0(algoParams);
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
%
%-------------------------
%
FF = outParams.FF;
fm = outParams.fm;
save(output_FF_file,'FF'); 
save(output_fm_file,'fm');
%
matrix_size = size(fm);
%-----------------------
if index_case <= 17 && algoParams.debug > 0
    load refdata.mat;
    ref_FF = REFCASES{index_case};
    clear MASKS;
    clear REFCASES;
    matrix_size = size(algoParams.complex_image);
    %size(ref_FF)
    %
    ref_3D = zeros(matrix_size(1:3));
    for index_slice = 1:matrix_size(3)
        for index_x = 1:matrix_size(2)
            for index_y = 1:matrix_size(1)
            jj = index_y + (index_x - 1)*matrix_size(1) + (index_slice - 1)*matrix_size(1)*matrix_size(2);
            ref_3D(index_y,index_x,index_slice) = ref_FF(jj,1);
            end
        end
    end
    %
%--------------------------------------------------------------------
    clear ratio
    for index_slice = 1:matrix_size(3)
    clear im_tmp;
    clear mask_slice;
    %
    mask_slice(:,:) = imDataParams.mask(:,:,index_slice);
    num_vec = mask_slice(:);
    num_vec(num_vec == 0) = [];
    im_tmp = outParams.FF(:,:,index_slice) - ref_3D(:,:,index_slice);
    im_tmp = im_tmp.*mask_slice;
    error_mask = zeros(matrix_size(1:2));
    error_mask(abs(im_tmp) >= 0.1) = 1;
    error_vec = error_mask(:);
    error_vec(error_vec == 0) = [];
    %
    ratio(index_slice) = 1- length(error_vec)/length(num_vec);
    end
    ratio
    mean(ratio)
    %
    for index_slice = 1:matrix_size(3)
    figure(1)
    subplot(1,matrix_size(3),index_slice);imagesc(fm(:,:,index_slice),[-400 400]);colormap gray;axis square;axis off;
    %
    figure(2)
    subplot(1,matrix_size(3),index_slice);imagesc(outParams.FF(:,:,index_slice),[0 1]);colormap gray;axis square;axis off;
    %
    figure(3)
    subplot(1,matrix_size(3),index_slice);imagesc(ref_3D(:,:,index_slice),[0 1]);colormap gray;axis square;axis off;    
    %
    %index_slice
    %pause
    end
    pause
elseif algoParams.debug > 0
    
    matrix_size = size(imDataParams.images);
%--------------------------------------------------------------------
    %
    for index_slice = 1:matrix_size(3)
    figure(1)
    subplot(1,matrix_size(3),index_slice);imagesc(fm(:,:,index_slice),[-400 400]);colormap gray;axis square;axis off;
    %
    figure(2)
    subplot(1,matrix_size(3),index_slice);imagesc(outParams.FF(:,:,index_slice),[0 1]);colormap gray;axis square;axis off;
    % 
    %
    %index_slice
    %pause
    end    
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%save(output_file,'outParams');
%write_raw(output_file,fm);
%-------------------------
%pause
index_case
%
end
