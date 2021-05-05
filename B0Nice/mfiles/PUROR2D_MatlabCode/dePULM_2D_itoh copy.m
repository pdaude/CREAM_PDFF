% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [phase_itoh_y mask_2D] = dePULM_2D_itoh(phase_tmp_ref,phase_tmp)
%----------------------------------------------------------------------
%   load the complex imaging file
    [yr_tmp xr_tmp] = size(phase_tmp);
%----------------------------------------------------------------------
%   Initializing phase_2D was done
%
phase_itoh_y = zeros(yr_tmp,xr_tmp);
mask_2D = zeros(yr_tmp,xr_tmp);
for index_x = 1:xr_tmp
        mag_index = find(phase_tmp_ref(:,index_x) ~= 0);
        th_index = find(phase_tmp_ref(:,index_x) == 0);
        mask_2D(mag_index,index_x) = 1;
        if (index_x/2 - round(index_x)/2)== 0
        phase_1D_tmp = phase_tmp(:,index_x);
        phase_itoh_y(:,index_x) = unwrap(phase_1D_tmp);
        else
        phase_1D_tmp = phase_tmp(xr_tmp:-1:1,index_x);
        phase_itoh_y(xr_tmp:-1:1,index_x) = unwrap(phase_1D_tmp);            
        end
        phase_itoh_y(th_index,index_x) = 0;
end