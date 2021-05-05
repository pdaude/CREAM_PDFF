%
function [algoParams] = CheckB0(algoParams)
pi_2 = 2*pi; 
%
complex_image = algoParams.complex_image;
%
TEratio = round(algoParams.delta_TE4B0/(algoParams.TE_seq(2)-algoParams.TE_seq(1)));
index_shift = -(TEratio-1):(TEratio-1);
if length(index_shift) == 1
index_shift = [-1 0 1];
end
index_shift(index_shift == 0) = [];
index_shift = [0 index_shift];
algoParams.index_shift = index_shift;
%-------------------------
matrix_size = algoParams.matrix_size;
%
%-------------------------
B0_map_ori = algoParams.B0map_Hz*(2*pi*algoParams.delta_TE4B0);
algoParams.B0map_Hz_ori = algoParams.B0map_Hz;
%
diffFF_ss = zeros(matrix_size(1), matrix_size(2),matrix_size(3),length(index_shift));
B0_map_shift = zeros(matrix_size(1), matrix_size(2),matrix_size(3),length(index_shift));
for index_ss = 1:length(index_shift)
    %-------------------
    B0_map_shift(:,:,:,index_ss) = B0_map_ori + index_shift(index_ss)*pi_2;
    %-------------------
    algoParams.B0map_Hz = algoParams.B0map_Hz_ori + index_shift(index_ss)/algoParams.delta_TE4B0;
    [algoParams] = FatCats_FFfromPhase(algoParams);
    %-------------------
    diffFF_ss(:,:,:,index_ss) =  (algoParams.phaseFF - algoParams.magFF);
    %-------------------
    %figure(1);imshow3D(algoParams.magFF,[-1 1]);colormap gray;
    %figure(2);imshow3D(algoParams.phaseFF,[-1 1]);colormap gray;
    %figure(3);imshow3D(diffFF_ss(:,:,:,index_ss),[-1 1]);colormap gray;
    %index_ss
    %pause
end
algoParams.B0map_Hz = algoParams.B0map_Hz_ori;
%
for index_slice = 1:matrix_size(3)
    BW_slice(:,:) = algoParams.BW_label(:,:,index_slice);
    %
    B0_out(:,:) = B0_map_ori(:,:,index_slice);
    %
    for index_BW = 1:max(BW_slice(:))
        for index_ss = 1:length(index_shift)
        diff_slice(:,:) = abs(diffFF_ss(:,:,index_slice,index_ss)).*abs(complex_image(:,:,index_slice,1,1));
        diff_slice(:,:) = diff_slice.*(1-algoParams.magFF(:,:,index_slice));
        diff_slice(isnan(diff_slice)) = 0;
            diff_slice(BW_slice ~= index_BW) = 0;
            sum_vec(index_ss) = sum(diff_slice(:));
        end
        %
        [mm Imin] = min(sum_vec);
        index_Imin = find(sum_vec == mm);
        II = index_Imin(1);
        if length(index_Imin) == 2 && abs(index_shift(index_Imin(1))) > abs(index_shift(index_Imin(2)))
        II = index_Imin(2);
        end
        %
        B0_ss(:,:) = B0_map_shift(:,:,index_slice,II);
        B0_out(BW_slice == index_BW) = B0_ss(BW_slice == index_BW);
        %
    end

    B0_map_ori(:,:,index_slice) = B0_out;

end

%
algoParams.B0map_Hz = B0_map_ori./(2*pi*algoParams.delta_TE4B0);
%-------------------------

