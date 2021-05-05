% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [phase_tmp] = dePULM_refY_trlX_tmp(phase_tmp_ref, phase_tmp,mask_mean,start_dw,start_up)
%----------------------------------------------------------------------
    pi_2 = pi*2.0;
    [yr_tmp xr_tmp] = size(phase_tmp);
%----------------------------------------------------------------------
    for index_y = 1:start_dw
        index_s = find(phase_tmp(index_y,:) ~= 0);
        index_s_mean = find(mask_mean(index_y,:) ~= 0);
        %--------------------
        if length(index_s_mean)/length(index_s) >= 0.3 && length(index_s_mean) > 6
        ref_tmp = phase_tmp_ref(index_y,index_s_mean);
        trl_tmp = phase_tmp(index_y,index_s_mean);
        else
        ref_tmp = phase_tmp_ref(index_y,index_s);
        trl_tmp = phase_tmp(index_y,index_s);            
        end
        %--------------------
        diff_tmp = mean_liu(ref_tmp - trl_tmp);
        %--------------------
        if abs(diff_tmp) > pi
        phase_tmp(index_y,index_s) = phase_tmp(index_y,index_s) + pi_2*round(diff_tmp/pi_2);     
        end
        %
    end   
    for index_y = start_up:yr_tmp
        index_s = find(phase_tmp(index_y,:) ~= 0);
        index_s_mean = find(mask_mean(index_y,:) ~= 0);
        %--------------------
        if length(index_s_mean)/length(index_s) >= 0.3 && length(index_s_mean) > 6
        ref_tmp = phase_tmp_ref(index_y,index_s_mean);
        trl_tmp = phase_tmp(index_y,index_s_mean);
        else
        ref_tmp = phase_tmp_ref(index_y,index_s);
        trl_tmp = phase_tmp(index_y,index_s);            
        end
        %--------------------
        diff_tmp = mean_liu(ref_tmp - trl_tmp);
        %--------------------
        if abs(diff_tmp) > pi
        phase_tmp(index_y,index_s) = phase_tmp(index_y,index_s) + pi_2*round(diff_tmp/pi_2);     
        end
        %
    end    
%---------------------------------------------------------------------- 
%   phase_unwrapping was done
