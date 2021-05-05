% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_2D] = dePULM_2D_xy_cen(phase_tmp_ref,phase_tmp,xy_start_dw,xy_start_up,mask_mean)
%----------------------------------------------------------------------
pi_2 = pi*2.0;
[yr_tmp xr_tmp] = size(phase_tmp);
%-----------------------------
    for index_x = 1:xr_tmp %index_x_min:index_x_max
        index_l = find(phase_tmp_ref(:,index_x) ~= 0);
        index_s_mean = find(mask_mean(xy_start_dw:xy_start_up,index_x) ~= 0);
        %--------------------
        if isempty(index_s_mean)
        index_s = find(mask_mean(:,index_x) ~= 0);            
        ref_tmp = phase_tmp_ref(index_s,index_x);
        trl_tmp = phase_tmp(index_s,index_x);
        %--------------------
        diff_tmp = mean_liu(ref_tmp - trl_tmp);
        shift_tmp = pi_2*round(diff_tmp/pi_2);
        phase_tmp(index_l,index_x) = phase_tmp(index_l,index_x) + shift_tmp;
        else
        index_s_mean = index_s_mean + xy_start_dw - 1;
        ref_tmp = phase_tmp_ref(index_s_mean,index_x);
        trl_tmp = phase_tmp(index_s_mean,index_x);
        %--------------------
        diff_tmp = mean_liu(ref_tmp - trl_tmp);
        %--------------------
            if abs(diff_tmp) > pi
            shift_tmp = pi_2*round(diff_tmp/pi_2);
            phase_tmp(index_l,index_x) = phase_tmp(index_l,index_x) + shift_tmp;
            end
        end
        %
    end 
%------------------------------------------------------------------------
out_2D = phase_tmp;%_xy;
%--------------------------------------------------------------------------
%   combination unwrapped_phase_x and _y was done
