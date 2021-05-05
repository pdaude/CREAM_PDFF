% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_2D] = dePULM_2D_diff_y(phase_tmp)
%----------------------------------------------------------------------
    pi_2 = pi*2.0;
    [yr_tmp xr_tmp] = size(phase_tmp);
%----------------------------------------------------------------------
    trl_dw = phase_tmp;
%------------
    for index_x = 2:(xr_tmp - 1)
        tmp_up = zeros(1,yr_tmp);
        tmp_dw = zeros(1,yr_tmp);
        tmp_tr = zeros(1,yr_tmp);
        tmp_up(trl_dw(:,index_x + 1) ~= 0) = 1;
        tmp_dw(trl_dw(:,index_x - 1) ~= 0) = 1;
        tmp_tr(trl_dw(:,index_x) ~= 0) = 1;
        index_x_b = find((tmp_up + tmp_dw + tmp_tr) == 3);            
        %---------------------
        diff_up = zeros(1,yr_tmp);
        diff_dw = zeros(1,yr_tmp);
        diff_up(index_x_b) =  trl_dw(index_x_b, index_x - 1) - trl_dw(index_x_b, index_x);
        diff_dw(index_x_b) =  trl_dw(index_x_b, index_x + 1) - trl_dw(index_x_b, index_x);
        q_up = find(abs(diff_up) > 1.0*pi);
        q_dw = find(abs(diff_dw) > 1.0*pi);
        if isempty(q_up) || isempty(q_dw)
        else
            quality_up_tmp = (length(index_x_b) - length(q_up))/length(index_x_b);
            quality_dw_tmp = (length(index_x_b) - length(q_dw))/length(index_x_b); 
            quality_br_tmp = quality_up_tmp*quality_dw_tmp; 
        %--------------------
            bad_up = zeros(1,yr_tmp);
            bad_dw = zeros(1,yr_tmp);
            if quality_br_tmp < 0.95
            %quality_up_tmp
            %quality_dw_tmp
         %index_y
            bad_up(q_up) = 2;
            bad_dw(q_dw) = 1;         
            bad_point = bad_up + bad_dw;
            fix_index = find(bad_point == 3);
            if isempty(fix_index)
            else
 %               tmp_up(index_x_b) =  trl_dw(index_x_b, index_x - 1) - trl_dw(index_x_b, index_x);
 %               tmp_dw(index_x_b) =  trl_dw(index_x_b, index_x + 1) - trl_dw(index_x_b, index_x);
                for ii = 1:length(fix_index)
                    jj = fix_index(ii);
                phase_tmp(jj, index_x) = phase_tmp(jj, index_x) +...
                                       pi_2*(round(0.5*(diff_up(jj)+ diff_dw(jj))/pi_2));
                end
            end
            end
        end
    end  
%---------------------------------------------------------------------
    index_x = 1;
        tmp_dw = zeros(1,yr_tmp);
        tmp_tr = zeros(1,yr_tmp);
        tmp_dw(phase_tmp(:, index_x + 1) ~= 0) = 1;
        tmp_tr(phase_tmp(:, index_x) ~= 0) = 1;
        index_y_b = find((tmp_dw + tmp_tr) == 2);            
        %-----------------------
        diff_dw = zeros(1,yr_tmp);
        diff_dw(index_y_b) = phase_tmp(index_y_b, index_x + 1) - phase_tmp(index_y_b, index_x);
        q_dw = find(abs(diff_dw) > 1.0*pi);
        if isempty(q_dw) || length(index_y_b) < 3
        else
            quality_dw_tmp = (length(index_y_b) - length(q_dw))/length(index_y_b); 
        %--------------------
            if quality_dw_tmp < 0.95
            fix_index = index_y_b;
                shift_tmp = [];
                shift_tmp(1:length(fix_index)) = phase_tmp(fix_index, index_x + 1) - phase_tmp(fix_index,index_x);
                phase_tmp(fix_index, index_x) = phase_tmp(fix_index, index_x) + pi_2*(round(shift_tmp'/pi_2));
            end
        end 
%---------------------------------------------------------------------
    index_x = xr_tmp;
        tmp_dw = zeros(1,yr_tmp);
        tmp_tr = zeros(1,yr_tmp);
        tmp_dw(phase_tmp(:,index_x - 1) ~= 0) = 1;
        tmp_tr(phase_tmp(:, index_x) ~= 0) = 1;
        index_y_b = find((tmp_dw + tmp_tr) == 2);            
        %-----------------------
        diff_dw = zeros(1,yr_tmp);
        diff_dw(index_y_b) = phase_tmp(index_y_b, index_x - 1) - phase_tmp(index_y_b, index_x);
        q_dw = find(abs(diff_dw) > 1.0*pi);
        if isempty(q_dw) || length(index_y_b) < 3
        else
            quality_dw_tmp = (length(index_y_b) - length(q_dw))/length(index_y_b); 
        %--------------------
            if quality_dw_tmp < 0.95
            fix_index = index_y_b;
                shift_tmp = [];
                shift_tmp(1:length(fix_index)) = phase_tmp(fix_index, index_x - 1) - phase_tmp(fix_index, index_x);
                phase_tmp(fix_index, index_x) = phase_tmp(fix_index, index_x) + pi_2*(round(shift_tmp'/pi_2));
            end
        end    
%---------------------------------------------------------------------
out_2D = phase_tmp;
%   phase_unwrapping was done
