% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_2D] = dePULM_2D_diff(phase_tmp)
%----------------------------------------------------------------------
    pi_2 = pi*2.0;
    [yr_tmp xr_tmp] = size(phase_tmp);
    %start_dw = round(yr_tmp/2);
    %start_up = start_dw + 1;
%----------------------------------------------------------------------
    trl_dw = phase_tmp;       
%------------
    for index_y = 2:(yr_tmp - 1)
        tmp_up = zeros(1,xr_tmp);
        tmp_dw = zeros(1,xr_tmp);
        tmp_tr = zeros(1,xr_tmp);
        tmp_up(trl_dw(index_y + 1,:) ~= 0) = 1;
        tmp_dw(trl_dw(index_y - 1,:) ~= 0) = 1;
        tmp_tr(trl_dw(index_y,:) ~= 0) = 1;
        index_x_b = find((tmp_up + tmp_dw + tmp_tr) == 3);            
        %-----------------------
        diff_up = zeros(1,xr_tmp);
        diff_dw = zeros(1,xr_tmp);
        diff_up(index_x_b) = trl_dw(index_y - 1,index_x_b) - trl_dw(index_y,index_x_b);
        diff_dw(index_x_b) = trl_dw(index_y + 1,index_x_b) - trl_dw(index_y,index_x_b);
        q_up = find(abs(diff_up) > 1.0*pi);
        q_dw = find(abs(diff_dw) > 1.0*pi);
        if isempty(q_up) || isempty(q_dw)
        else
            quality_up_tmp = (length(index_x_b) - length(q_up))/length(index_x_b);
            quality_dw_tmp = (length(index_x_b) - length(q_dw))/length(index_x_b); 
            quality_br_tmp = quality_up_tmp*quality_dw_tmp; 
        %--------------------
            bad_up = zeros(1,xr_tmp);
            bad_dw = zeros(1,xr_tmp);
            if quality_br_tmp < 0.95
            bad_up(q_up) = 2;
            bad_dw(q_dw) = 1;         
            bad_point = bad_up + bad_dw;
            fix_index = find(bad_point == 3);
                if isempty(fix_index)
                else
                shift_tmp = [];
                shift_tmp(1:length(fix_index)) = (trl_dw(index_y + 1,fix_index) + trl_dw(index_y - 1,fix_index)...
                                               - 2.0*trl_dw(index_y,fix_index))/2;
                phase_tmp(index_y,fix_index) = phase_tmp(index_y,fix_index) + pi_2*(round(shift_tmp/pi_2));
                end
            end
        end
    end   
%---------------------------------------------------------------------
    index_y = 1;
        tmp_dw = zeros(1,xr_tmp);
        tmp_tr = zeros(1,xr_tmp);
        tmp_dw(phase_tmp(index_y + 1,:) ~= 0) = 1;
        tmp_tr(phase_tmp(index_y,:) ~= 0) = 1;
        index_x_b = find((tmp_dw + tmp_tr) == 2);            
        %-----------------------
        diff_dw = zeros(1,xr_tmp);
        diff_dw(index_x_b) = phase_tmp(index_y + 1,index_x_b) - phase_tmp(index_y,index_x_b);
        q_dw = find(abs(diff_dw) > 1.0*pi);
        if isempty(q_dw) || length(index_x_b) < 3
        else
            quality_dw_tmp = (length(index_x_b) - length(q_dw))/length(index_x_b); 
        %--------------------
            if quality_dw_tmp < 0.95
            fix_index = index_x_b;
                shift_tmp = [];
                shift_tmp(1:length(fix_index)) = phase_tmp(index_y + 1,fix_index) - phase_tmp(index_y,fix_index);
                phase_tmp(index_y,fix_index) = phase_tmp(index_y,fix_index) + pi_2*(round(shift_tmp/pi_2));
            end
        end 
%---------------------------------------------------------------------
    index_y = yr_tmp;
        tmp_dw = zeros(1,xr_tmp);
        tmp_tr = zeros(1,xr_tmp);
        tmp_dw(phase_tmp(index_y - 1,:) ~= 0) = 1;
        tmp_tr(phase_tmp(index_y,:) ~= 0) = 1;
        index_x_b = find((tmp_dw + tmp_tr) == 2);            
        %-----------------------
        diff_dw = zeros(1,xr_tmp);
        diff_dw(index_x_b) = phase_tmp(index_y - 1,index_x_b) - phase_tmp(index_y,index_x_b);
        q_dw = find(abs(diff_dw) > 1.0*pi);
        if isempty(q_dw) || length(index_x_b) < 3
        else
            quality_dw_tmp = (length(index_x_b) - length(q_dw))/length(index_x_b); 
        %--------------------
            if quality_dw_tmp < 0.95
            fix_index = index_x_b;
                shift_tmp = [];
                shift_tmp(1:length(fix_index)) = phase_tmp(index_y - 1,fix_index) - phase_tmp(index_y,fix_index);
                phase_tmp(index_y,fix_index) = phase_tmp(index_y,fix_index) + pi_2*(round(shift_tmp/pi_2));
            end
        end    
%---------------------------------------------------------------------
out_2D = phase_tmp;
%   phase_unwrapping was done
