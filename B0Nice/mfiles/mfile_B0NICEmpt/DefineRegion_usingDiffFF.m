%
function [algoParams] = DefineRegion_usingDiffFF(algoParams)
%
se = strel('disk',2);
%
[algoParams] = FatCats_FFfromPhase(algoParams);
%
Res_map = algoParams.phaseFF;
matrix_size = algoParams.matrix_size;
BW_label = zeros(matrix_size(1:3));
for index_slice = 1:matrix_size(3)
    %
    Res_tmp(:,:) = Res_map(:,:,index_slice);
    %----------------------------
    Res_tmp(abs(algoParams.phaseFF(:,:,index_slice) - algoParams.magFF(:,:,index_slice)) < 0.75) = 0;
    Res_tmp(Res_tmp ~= 0) = 1;
    %
    [L, num] = bwlabel(Res_tmp, 4); 
    s  = regionprops(L, 'area');
    [y_area I_area] = sort([s.Area],'descend');

    y_area(y_area <algoParams.minAreaUsingFlip) = [];

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
    %subplot(1,2,1);imagesc(Res_tmp);colormap gray;
    %subplot(1,2,2);imagesc(BW_label(:,:,index_slice),[0 6]);colormap gray;
    %index_slice
    %pause
end
algoParams.BW_label = BW_label;




