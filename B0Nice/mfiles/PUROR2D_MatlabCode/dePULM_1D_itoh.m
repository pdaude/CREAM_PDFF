% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [phase_itoh_x phase_itoh_y] = dePULM_1D_itoh(phase_original)
%----------------------------------------------------------------------
%   load the complex imaging file
    [yr_tmp xr_tmp] = size(phase_original);
%----------------------------------------------------------------------

%   Initializing phase_2D was done
    phase_itoh_x = zeros(yr_tmp,xr_tmp);
for index_y = 1:yr_tmp
    if mode(index_y,2) == 0
    phase_1D_tmp = phase_original(index_y,:);
    phase_itoh_x(index_y,:) = unwrap(phase_1D_tmp);
    else
    phase_1D_tmp = phase_original(index_y,xr_tmp:-1:1);
    phase_itoh_x(index_y,xr_tmp:-1:1) = unwrap(phase_1D_tmp);        
    end
end
%
phase_itoh_y = zeros(yr_tmp,xr_tmp);
for index_x = 1:xr_tmp
        if mode(index_x,2) == 0
        phase_1D_tmp = phase_original(:,index_x);
        phase_itoh_y(:,index_x) = unwrap(phase_1D_tmp);
        else
        phase_1D_tmp = phase_original(yr_tmp:-1:1,index_x);
        phase_itoh_y(yr_tmp:-1:1,index_x) = unwrap(phase_1D_tmp);            
        end
end