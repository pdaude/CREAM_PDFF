% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [quality_br xy_start xy_end] = dePULM_2D_quality_y(phase_tmp,mask_q)
%----------------------------------------------------------------------
%    
%   check the quality of the bridge pixels
   [yr_tmp xr_tmp] = size(phase_tmp);
   quality_up = zeros(1,yr_tmp);
   quality_dw = zeros(1,yr_tmp);
   quality_br = -ones(1,yr_tmp);
    for index_x = 1:xr_tmp
        if index_x ~= 1 && index_x ~= xr_tmp
            tmp_up = zeros(1,yr_tmp);
            tmp_dw = zeros(1,yr_tmp);
            tmp_tr = zeros(1,yr_tmp);
            tmp_up(mask_q(:,index_x + 1) ~= 0) = 1;
            tmp_dw(mask_q(:,index_x - 1) ~= 0) = 1;
            tmp_tr(mask_q(:,index_x) ~= 0) = 1;
            index_q = find((tmp_up + tmp_dw + tmp_tr) == 3);
            if isempty(index_q)
            else
            len_q = length(index_q);
            diff_up = phase_tmp(index_q,index_x) -  phase_tmp(index_q, index_x - 1);
            diff_dw = phase_tmp(index_q,index_x) -  phase_tmp(index_q, index_x + 1);
                q_up = find(abs(diff_up) > 1.0*pi);
                if isempty(q_up)
                len_up = 0;
                else
                len_up = length(q_up);
                end
        %    
                q_dw = find(abs(diff_dw) > 1.0*pi);
                if isempty(q_dw)
                len_dw = 0;
                else
            	len_dw = length(q_dw);
                end
            quality_up(index_x) = (len_q - len_up)/len_q;
            quality_dw(index_x) = (len_q - len_dw)/len_q;
            quality_br(index_x) = quality_up(index_x)*quality_dw(index_x);
            end
        else
            if index_x == 1
            tmp_up = zeros(1,yr_tmp);
            tmp_tr = zeros(1,yr_tmp);
            tmp_up(mask_q(:,index_x + 1) ~= 0) = 1;
            tmp_tr(mask_q(:,index_x) ~= 0) = 1;
            index_q = find((tmp_up + tmp_tr) == 2);
                if isempty(index_x)
                else
                len_q = length(index_q);
                diff_up = phase_tmp(index_q,index_x) -  phase_tmp(index_q, index_x + 1);
                q_up = find(abs(diff_up) > 1.0*pi);
                if isempty(q_up)
                len_up = 0;
                else
                len_up = length(q_up);
                end
                %    
                quality_br(index_x) = (len_q - len_up)/len_q;
                end            
            else
            tmp_up = zeros(1,yr_tmp);
            tmp_tr = zeros(1,yr_tmp);
            tmp_up(mask_q(:,index_x - 1) ~= 0) = 1;
            tmp_tr(mask_q(:,index_x) ~= 0) = 1;
            index_q = find((tmp_up + tmp_tr) == 2);
                if isempty(index_x)
                else
                len_q = length(index_q);
                diff_up = phase_tmp(index_q,index_x) -  phase_tmp(index_q,index_x - 1);
                q_up = find(abs(diff_up) > 1.0*pi);
                if isempty(q_up)
                len_up = 0;
                else
                len_up = length(q_up);
                end
                %    
                quality_br(index_x) = (len_q - len_up)/len_q;
                end                 
            end
        end
    end
%   
%    figure(8)
%    plot(quality_up,'o-b');hold on;
%    plot(quality_dw,'o-r');hold on;
%    plot(quality_br,'o-g');hold off;
%    pause
    %------------------------------------------------------------------------
%   find the protect region
    xy_start = 1;
    xy_end = 1;
    flag_start = 1;
    flag_end = 2;
    q_th = 0.9;
    for ss = 1:xr_tmp
        if quality_br(ss) > q_th && flag_start == 1
            g_start = ss;
            flag_start = 0;
            flag_end = 1;
        end
        %
        if quality_br(ss)< q_th && flag_end == 1
            g_end = ss - 1;
            flag_start = 1;
            flag_end = 0;
        end
        %
        if ss == xr_tmp && flag_end == 1
            if quality_br(ss) >= q_th
            g_end = ss;
            flag_end = 0;
            else
            g_end = g_start;
            flag_end = 0;                
            end
        end
        %
        if flag_end == 0 && ss > 1
           if (g_end - g_start) >= (xy_end - xy_start)
           xy_start = g_start;
           xy_end = g_end;
           end
        end
        %
    end
   %-----------------

