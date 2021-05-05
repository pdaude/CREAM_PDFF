%
function [algoParams] = DefineRegion_usingPhaseGradient(algoParams)
%
%
B0_map = (2*pi*algoParams.delta_TE4B0).*algoParams.B0map_Hz;
%-------------------------
th4label = algoParams.th4RegionDivision;%
BW_label = zeros(algoParams.matrix_size(1:3));
%------------------------
    se = strel('disk',2);        
for index_slice = 1:algoParams.matrix_size(3)
    %---------------------------
    mask_x = zeros(algoParams.matrix_size(1:2));
    mask_y = zeros(algoParams.matrix_size(1:2));
    %---------------------------
    B0_unit(:,:) = B0_map(:,:,index_slice);
    [Gx, Gy] = imgradientxy(B0_unit,'central');
    mask_x(abs(Gx) < th4label) = 1;
    mask_x = imerode(mask_x,se);
    %----------------------------    
    mask_y(abs(Gy) < th4label) = 1;
    mask_y = imerode(mask_y,se);
    %----------------------------
    mask_xy = mask_x.*mask_y;
    %
    [L, num] = bwlabel(mask_xy, 4);   

    s  = regionprops(L, 'area');
    [y_area I_area] = sort([s.Area],'descend');

    y_area(y_area <algoParams.minAreaUsingPhaseGradient) = [];

    BW_L_tmp = zeros(algoParams.matrix_size(1:2));
    for index_L = 1:length(y_area)
        BW_tmp = zeros(algoParams.matrix_size(1:2));
        BW_tmp(L == I_area(index_L)) = 1;
        %
        BW_tmp = imdilate(BW_tmp,se);
        BW_L_tmp(BW_tmp == 1) = index_L;
    end
    %
    BW_L_tmp = imdilate(BW_L_tmp,se);
    BW_label(:,:,index_slice) = BW_L_tmp;
    %
    if algoParams.flag_debug == 1
    figure(8)
    imagesc(BW_label(:,:,index_slice),[0 5]);axis square;axis off;
    figure(9)
    %subplot(1,3,2);imagesc(algoParams.FF_bin(:,:,index_slice));axis square;axis off; 
    imagesc(algoParams.B0_map(:,:,index_slice));axis square;axis off;     
    index_slice
    pause
    end
    %
end
algoParams.BW_label = BW_label;
%

%
