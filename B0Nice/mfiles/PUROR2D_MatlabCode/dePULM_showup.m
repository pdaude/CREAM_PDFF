% This program is to show the final results of the unwrapped phase
function [out_phase] = dePULM_showup(complex_image, line_unwrapping,slice_start,slice_end,line_mask)
%
%--------------------------------------------
    for index_slice = slice_start:slice_end
        index_slice
        %
        TE = 0.0217;
        size_filter = 13;
        mask = line_mask{index_slice};
        phase_tmp = line_unwrapping{index_slice}(:,:);
        [yr_tmp xr_tmp] = size(phase_tmp);
        figure(45)
        subplot(2,2,2);
        IM_mask = zeros(yr_tmp,xr_tmp);
        IM_mask(mask == 1) = phase_tmp(mask == 1);
        out_phase{index_slice} = IM_mask;
        imshow(IM_mask,[-4.0*pi  4.0*pi]);
        % 
    %---------------------------------------
    %
    unwrapped_phase = phase_tmp;
    H = fspecial('average',size_filter);
    IM = unwrapped_phase - imfilter(unwrapped_phase,H);
    IM_mask = zeros(yr_tmp,xr_tmp);
    IM_mask(mask == 1) = IM(mask == 1);
    IM_mask = IM_mask/(2.0*pi*TE);
    subplot(2,2,4)
    imshow(-IM_mask,[-6 6]);
    david_combined_2D_filtered{index_slice} = IM_mask;
    %tmp = unwrapped_phase(:);
    %tmp = angle(exp(i.*tmp));
    %unwrapped_phase = reshape(tmp,yr_tmp,xr_tmp);
    %------------------------------------------
    %imshow(unwrapped_phase,[-0.22*pi 0.22*pi]);
    %
    size_tmp = size(complex_image);
    if length(size_tmp) < 3
    mag_original = abs(complex_image{index_slice}(:,:));
    phase_original = angle(complex_image{index_slice}(:,:));
    else
    mag_original = abs(complex_image(:,:,index_slice));
    phase_original = angle(complex_image(:,:,index_slice));
    end
    subplot(2,2,3)
    imagesc(mag_original);colormap gray;axis square;axis off;
    %
    subplot(2,2,1)
    imshow(phase_original,[-1.0*pi  1.0*pi]);
    %--------------------------------------
    pause
    end
    unwrapped_phase_3D = cat(3, line_unwrapping{:});
    H = fspecial3('gaussian',[size_filter size_filter 3]);
    D3_new = imfilter(unwrapped_phase_3D,H,'replicate');
    size(D3_new)
    for index_slice = slice_start:slice_end
    mag_original = abs(complex_image{index_slice}(:,:));
    phase_original = angle(complex_image{index_slice}(:,:));
    subplot(2,2,1)
    imshow(phase_original,[-1.0*pi  1.0*pi]);   
    subplot(2,2,3)
    imagesc(mag_original);colormap gray;axis square;axis off;
    subplot(2,2,2)
    imshow(unwrapped_phase_3D(:,:,index_slice),[-4.0*pi  4.0*pi]);    
    %---------------------------------------------------------
        IM = unwrapped_phase_3D(:,:,index_slice) - D3_new(:,:,index_slice);
        mask = line_mask{index_slice};
        IM_mask = zeros(yr_tmp,xr_tmp);
        IM_mask(mask == 1) = IM(mask == 1);
        IM_mask = IM_mask/(2.0*pi*TE);
        david_combined_3D_filtered{index_slice} = IM_mask;
        figure(45)
        subplot(2,2,4)        
    imshow(-IM_mask,[-6 6]);
    pause
    end 