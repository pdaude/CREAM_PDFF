% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_2D] = dePULM_2D_xy(unwrapped_phase_x,unwrapped_phase_y,xy_start_dw,xy_start_up,index_min_max)
%----------------------------------------------------------------------
index_y_min = index_min_max(1);
index_y_max = index_min_max(2);
index_x_min = index_min_max(3);
index_x_max = index_min_max(4);
phase_tmp_xy = unwrapped_phase_x;
%-----------------------------
st_phase_x = unwrapped_phase_x(xy_start_dw,:) - unwrapped_phase_y(xy_start_dw,:);
for index_y = xy_start_dw:-1:index_y_min
st_tmp = zeros(1,length(st_phase_x));
st_tmp(unwrapped_phase_x(index_y,:) ~= 0) = 1;
phase_tmp_xy(index_y,index_x_min:index_x_max) = unwrapped_phase_y(index_y,index_x_min:index_x_max)...
                                           + st_phase_x(index_x_min:index_x_max).*st_tmp(index_x_min:index_x_max);
end                                                            
%
st_phase_x = unwrapped_phase_x(xy_start_up,:) - unwrapped_phase_y(xy_start_up,:);
for index_y = xy_start_up:index_y_max
st_tmp = zeros(1,length(st_phase_x));
st_tmp(unwrapped_phase_x(index_y,:) ~= 0) = 1;
phase_tmp_xy(index_y,index_x_min:index_x_max) = unwrapped_phase_y(index_y,index_x_min:index_x_max)...
                                           + st_phase_x(index_x_min:index_x_max).*st_tmp(index_x_min:index_x_max);
end
%------------------------------------------------------------------------
out_2D = phase_tmp_xy;
%--------------------------------------------------------------------------
%   combination unwrapped_phase_x and _y was done
