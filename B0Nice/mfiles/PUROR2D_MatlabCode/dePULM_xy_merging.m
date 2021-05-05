% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_xy] = dePULM_xy_merging(out_phase_x, out_phase_y,cell_signal, cell_connect, index_min_max,start_dw,start_up)
%----------------------------------------------------------------------
%   load the imaging file
    index_y_min = index_min_max(1);
    index_y_max = index_min_max(2);
%----------------   
    pi_2 = pi*2;
    phase_range = 4.0*pi; %[-phase_range phase_range]
    [yr xr] = size(out_phase_x);
    jump_w = pi;
    out_xy = out_phase_x;
%--------------------------------------------------------------
    q_out_xy = zeros(1,yr);
    q_out_xy(start_dw:start_up) = 1;
    %
    for index_y = start_dw:-1:index_y_min
        index_s = cell_signal{1,index_y};
        diff_tmp = out_phase_y(index_y,index_s) - out_phase_x(index_y,index_s);
        index_j = diff_tmp(abs(diff_tmp) <= pi);
        if isempty(index_j)
            q_out_tmp = 0;
        else
        q_out_tmp = length(index_j)/length(index_s);
        end
        diff_q = q_out_xy(index_y + 1) - q_out_tmp;
        if diff_q > 0.2
           phase_RF = out_phase_y(index_y,index_s);
           phase_TR = out_phase_x(index_y,index_s);
           [index_ls_RF] = dePULM_2D_ls(phase_RF,index_s);
        %    
           ls_mag = index_s;
           index_ls = index_ls_RF;
           line_phase = phase_TR;        
           re_phase = phase_RF;
           [out_TR] = dePULM_2D_merging(line_phase,index_ls,re_phase);
           out_xy(index_y,ls_mag) = out_TR;
           %
           % if q_out_tmp <0.5
           % out_xy(index_y,:) = out_phase_y(index_y,:);
           % q_out_xy(index_y) = 1;
           % else
           % out_xy(index_y,index_s) = out_xy(index_y,index_s) + pi_2*round(diff_tmp/(pi_2*1.1));    
            diff_tmp = out_phase_y(index_y,index_s) - out_xy(index_y,index_s);
            index_j = diff_tmp(abs(diff_tmp) <= pi);
            q_out_xy(index_y) = length(index_j)/length(index_s);           
            %end
        else
            q_out_xy(index_y) = q_out_tmp;
        end
        
    end
    %
    for index_y = start_up:index_y_max
        index_s = cell_signal{1,index_y};
        diff_tmp = out_phase_y(index_y,index_s) - out_phase_x(index_y,index_s);
        index_j = diff_tmp(abs(diff_tmp) <= pi);
        if isempty(index_j)
            q_out_tmp = 0;
        else
        q_out_tmp = length(index_j)/length(index_s);
        end
        %
        diff_q = q_out_xy(index_y - 1) - q_out_tmp;
        if diff_q > 0.2
           phase_RF = out_phase_y(index_y,index_s);
           phase_TR = out_phase_x(index_y,index_s);
           [index_ls_RF] = dePULM_2D_ls(phase_RF,index_s); 
        %    
           ls_mag = index_s;
           index_ls = index_ls_RF;
           line_phase = phase_TR;        
           re_phase = phase_RF;
           [out_TR] = dePULM_2D_merging(line_phase,index_ls,re_phase);
           out_xy(index_y,ls_mag) = out_TR;            
           diff_tmp = out_phase_y(index_y,index_s) - out_xy(index_y,index_s);
           index_j = diff_tmp(abs(diff_tmp) <= pi);
           q_out_xy(index_y) = length(index_j)/length(index_s);
        else
            q_out_xy(index_y) = q_out_tmp;
        end     
    end    
%    figure(3)
%    plot(q_out_xy,'o-r');
%---------------------------------------------------------
%--------------------------------------------------------------------------
%   Combination of out_phase_x and _y was done
