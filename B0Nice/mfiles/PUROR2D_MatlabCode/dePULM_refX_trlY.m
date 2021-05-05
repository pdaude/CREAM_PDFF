% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_2D] = dePULM_refX_trlY(phase_tmp_ref, phase_tmp,index_min_max, start_dw,start_up)
%----------------------------------------------------------------------
    pi_2 = pi*2.0;
    index_y_min = index_min_max(1);
    index_y_max = index_min_max(2);
    index_x_min = index_min_max(3);
    index_x_max = index_min_max(4);
%----------------------------------------------------------------------
    ref_dw = phase_tmp_ref(index_y_min:start_dw,:);
    trl_dw = phase_tmp(index_y_min:start_dw,:);
    %--------
    ref_up = phase_tmp_ref(start_up:index_y_max,:);
    trl_up = phase_tmp(start_up:index_y_max,:);
%------------
    for index_x = index_x_min:index_x_max
        index_s_dw = find(ref_dw(:,index_x) ~= 0);
        index_s_up = find(ref_up(:,index_x) ~= 0);
        %--------------------
        ref_dw_tmp = ref_dw(index_s_dw,index_x);
        trl_dw_tmp = trl_dw(index_s_dw,index_x);
        %--------------------
        diff_tmp = mean_liu(ref_dw_tmp - trl_dw_tmp);
        %--------------------
        if abs(diff_tmp) > pi
        trl_dw(index_s_dw,index_x) = trl_dw_tmp + pi_2*round(diff_tmp/pi_2);
        trl_dw_tmp = trl_dw(index_s_dw,index_x);
        end
        %   
        out_TR = trl_dw_tmp;
        %
        trl_dw(index_s_dw,index_x) = out_TR;
%--------------
        ref_up_tmp = ref_up(index_s_up,index_x);
        trl_up_tmp = trl_up(index_s_up,index_x);
        %--------------------
        diff_tmp = mean_liu(ref_up_tmp - trl_up_tmp);
        %--------------------        
        if abs(diff_tmp) > pi
        trl_up(index_s_up,index_x) = trl_up_tmp + pi_2*round(diff_tmp/pi_2);
        trl_up_tmp = trl_up(index_s_up,index_x);
        end
        %--------------------   
        out_TR = trl_up_tmp;
        %
        trl_up(index_s_up,index_x) = out_TR;        
    end    
%---------------------------------------------------------------------- 
    phase_tmp(index_y_min:start_dw,:) = trl_dw;
    phase_tmp(start_up:index_y_max,:) = trl_up;
%-----------------
out_2D = phase_tmp;
%   phase_unwrapping was done
