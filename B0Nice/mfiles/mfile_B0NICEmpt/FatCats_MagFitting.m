function [algoParams] = FatCats_MagFitting(algoParams)
%
FF_trl = algoParams.FF_trl;
%
    debug = algoParams.debug;
%    
    matrix_size = size(algoParams.complex_image);
    if matrix_size(4) > 1
        mag = abs(algoParams.complex_image).^2;
        sum_mag = sqrt(sum(mag,4));
        algoParams.mag4magfitting = sum_mag;
    else
        algoParams.mag4magfitting = abs(algoParams.complex_image);        
    end
    algoParams.TE4magfitting = algoParams.TE_seq;
    %
    if length(algoParams.TE_seq) > algoParams.CutOff
        algoParams.mag4magfitting(:,:,:,:,(algoParams.CutOff+1):length(algoParams.TE_seq)) = [];
        algoParams.TE4magfitting((algoParams.CutOff+1):length(algoParams.TE_seq)) = [];
    end
    %
complex_image = algoParams.mag4magfitting;
TE_seq = algoParams.TE4magfitting;
%
model_r = algoParams.model_r;
model_f = algoParams.model_f;
index_min = 1;%algoParams.index_min;
matrix_size = size(complex_image);
%
    for index_trl = 1:length(FF_trl)
        FF_R2star = FF_trl(index_trl)*ones(matrix_size(1),matrix_size(2),matrix_size(3));
        [R2star_map_tmp R2star_residue_map] = FatCats_T2starMap(model_r,model_f,TE_seq,index_min,complex_image,matrix_size,FF_R2star,algoParams.R2star_pt,debug);
        R2star_map_tmp(isnan(R2star_map_tmp)) = 0;
        R2star_residue_map(isnan(R2star_residue_map)) = 1;
        if index_trl == 1
        R2star_weight = 1./(R2star_residue_map);
        R2star_map_raw = R2star_map_tmp./R2star_residue_map;
        FF_raw_mean = FF_trl(index_trl)./(R2star_residue_map);
        %
        else
        R2star_map_raw = R2star_map_raw + R2star_map_tmp./R2star_residue_map;
        R2star_weight = R2star_weight + 1./(R2star_residue_map);
        FF_raw_mean = FF_raw_mean + FF_trl(index_trl)./(R2star_residue_map);
        end    
    %
    %index_trl
    res_FF(:,:,:,index_trl) = R2star_residue_map;
    end
%
res_ratio(:,:,:) = (res_FF(:,:,:,1))./(res_FF(:,:,:,length(FF_trl)));
%
R2star_map_raw = R2star_map_raw./R2star_weight;
R2star_map_raw(isnan(R2star_map_raw)) = 0;
algoParams.R2star_map_raw = R2star_map_raw;
%
FF_raw_mean = FF_raw_mean./R2star_weight;
FF_raw_mean(isnan(FF_raw_mean)) = 0.5;
algoParams.magFF = FF_raw_mean;
%
    if debug == 1
    for index_slice = 1:matrix_size(3)
    im_tmp(:,:) = FF_raw_mean(:,:,index_slice);
    figure(1)
    subplot(1,matrix_size(3),index_slice);imagesc(im_tmp,[0 1]);axis square; axis off;
    %
    figure(2)
    subplot(1,matrix_size(3),index_slice);imagesc(R2star_map_raw(:,:,index_slice),[0 300]);colormap gray;axis square; axis off;
    %
    figure(3)
    subplot(1,matrix_size(3),index_slice);imagesc(res_ratio(:,:,index_slice),[0 2]);colormap gray;axis square; axis off;
    %index_slice
    %pause
    end
    pause
    end
%
end
%
function [R2star_map R2star_residue_map] = FatCats_T2starMap(model_r,model_f,TE_seq,index_min,complex_image,matrix_size,FF_map,R2star_pt,debug)
%
%--------------------------
index_start = index_min;
pi_2 =2*pi;
%---------------------------------------------------------
S_10 = abs(complex_image);
%-------------------------------------------------------
TE_0 = TE_seq(index_start);
C_TE_0 = model_r.*exp(1i*(2*pi*TE_0.*model_f));
C_tmp_0 = sum(C_TE_0);
for index_echo = 1:length(TE_seq)
    TE_1 = TE_seq(index_echo);
    C_TE_1 = model_r.*exp(1i*(2*pi*TE_1.*model_f));
    C_tmp_1(index_echo) = sum(C_TE_1);    
end
%
R2star_map = zeros(matrix_size(1),matrix_size(2),matrix_size(3));
R2star_residue_map = zeros(matrix_size(1),matrix_size(2),matrix_size(3));
for index_slice = 1:matrix_size(3)
    for index_y = 1:matrix_size(1)
        %for index_x = 1:matrix_size(2)
            S10_pt = S_10(index_y,:,index_slice,:,:);
            S10_pt = squeeze(S10_pt);
            %
            %B0_pt = B0_map(index_y,index_x,index_slice);
            %
            FF_percent = FF_map(index_y,:,index_slice);
            FF_percent = squeeze(FF_percent);
            [matrix_FF matrix_R2] = meshgrid(FF_percent,R2star_pt);
            %
            %calculate the R2star map
            residue_R2star = zeros(size(matrix_FF));
            S10_line = zeros(size(matrix_FF));
            for index_echo = 1:length(TE_seq)% (index_start+1):length(TE_seq)
            R2star_tmp = exp(-matrix_R2*(TE_seq(index_echo) - TE_seq(index_start))); 
            %
            FFeffect = abs((1-matrix_FF) + matrix_FF*C_tmp_1(index_echo))./abs((1-matrix_FF) + matrix_FF*C_tmp_0);
            %
                S10_pt_tmp = S10_pt(:,index_echo)./S10_pt(:,index_start);
                S10_pt_tmp(isnan(S10_pt_tmp)) = 1;
                %
                for index_R2 = 1:length(R2star_pt)
                S10_line(index_R2,:) = squeeze(S10_pt_tmp);  
                end
            %
                if index_echo ~= index_start %|| abs(phase_echo(index_echo)) < 0.9*pi
                residue_R2star = residue_R2star + ((FFeffect.*R2star_tmp - S10_line).^2);
                end
            end
            %
            [residue_R2star_tmp min_index_tmp] = min(residue_R2star,[],1);
            %
            R2star_map(index_y,:,index_slice) = R2star_pt(min_index_tmp);     
            R2star_residue_map(index_y,:,index_slice) = residue_R2star_tmp;
            %
        %end %index_x
    end
    R2star_slice(:,:) = R2star_map(:,:,index_slice);
    R2star_slice = medfilt2(R2star_slice,[3 3],'symmetric');
    R2star_map(:,:,index_slice) = R2star_slice;
    %close all
    %imagesc(R2star_slice,[0 300]);colormap gray;
    %index_slice
    %pause
end
%-------------------------------------------------------
if debug == 1
for index_slice = 1:matrix_size(3)
    R2star_slice(:,:) = R2star_map(:,:,index_slice);
    R2star_residue_slice(:,:) = R2star_residue_map(:,:,index_slice);
    %
    figure(12)
    subplot(1,2,1);imagesc(R2star_slice,[0 100]);colormap jet;axis square;axis off;
    subplot(1,2,2);imagesc(R2star_residue_slice,[0 0.5]);colormap jet;axis square;axis off;
    %
    index_slice
    pause
end
end
%------------------------------------------------------
end