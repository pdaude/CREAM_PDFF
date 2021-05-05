% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_2D] = dePULM_refX_trlY_tmp(phase_tmp_ref, phase_tmp,index_min_max,mask_mean,xy_start_dw,xy_start_up)
%----------------------------------------------------------------------
    pi_2 = pi*2.0;
    index_y_min = index_min_max(1);
    index_y_max = index_min_max(2);
    index_x_min = index_min_max(3);
    index_x_max = index_min_max(4);
%----------------------------------------------------------------------
    for index_x = index_x_min:index_x_max
        index_s = find(phase_tmp_ref(index_y_min:xy_start_dw,index_x) ~= 0);
        index_s_mean = find(mask_mean(index_y_min:xy_start_dw,index_x) ~= 0);
        %--------------------
        ref_tmp = phase_tmp_ref(index_s_mean,index_x);
        trl_tmp = phase_tmp(index_s_mean,index_x);
        %--------------------
        diff_tmp = mean_liu(ref_tmp - trl_tmp);
        %--------------------
        if abs(diff_tmp) > pi
        phase_tmp(index_s,index_x) = phase_tmp(index_s,index_x) + pi_2*round(diff_tmp/pi_2);     
        end
        %------------------
        index_s = find(phase_tmp_ref(xy_start_up:index_y_max,index_x) ~= 0);
        index_s_mean = find(mask_mean(xy_start_up:index_y_max,index_x) ~= 0);        
        %--------------------
        ref_tmp = phase_tmp_ref(index_s_mean,index_x);
        trl_tmp = phase_tmp(index_s_mean,index_x);
        %--------------------
        diff_tmp = mean_liu(ref_tmp - trl_tmp);
        %--------------------
        if abs(diff_tmp) > pi
        phase_tmp(index_s,index_x) = phase_tmp(index_s,index_x) + pi_2*round(diff_tmp/pi_2);   
        end
        %       
    end    
%---------------------------------------------------------------------- 
%-----------------
out_2D = phase_tmp;
%   phase_unwrapping was done
