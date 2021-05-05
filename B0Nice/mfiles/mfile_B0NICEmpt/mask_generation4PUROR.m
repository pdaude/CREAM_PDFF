%
function [algoParams] = mask_generation4PUROR(algoParams)
%-------------------------
matrix_size = algoParams.matrix_size;
%
   algoParams.mask4unwrap = zeros(matrix_size(1),matrix_size(2),matrix_size(3));
   algoParams.mask4supp = algoParams.mask4unwrap;
   algoParams.mask4STAS = algoParams.mask4unwrap;
   for index_slice = 1:matrix_size(3)
   mag_slice(:,:) = abs(algoParams.complex_B0(:,:,index_slice));     
   %--------------------------------------
        th_tmp = prctile(mag_slice(:),[algoParams.th4unwrap algoParams.th4supp algoParams.th4STAS]);
        algoParams.mask4unwrap(:,:,index_slice) = mag_slice > th_tmp(1);
   %---------------------------------------     
        algoParams.mask4supp(:,:,index_slice) = mag_slice > th_tmp(2);
   %---------------------------------------     
        algoParams.mask4STAS(:,:,index_slice) = mag_slice > th_tmp(3);
   %---------------------------------------  
   end
%
