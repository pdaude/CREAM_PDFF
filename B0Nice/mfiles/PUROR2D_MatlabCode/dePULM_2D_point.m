% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_2D] = dePULM_2D_point(phase_tmp_ref, phase_tmp,index_min_max, start_dw,start_up,xr)
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
    for index_y = (start_dw - index_y_min - 0):-1:2
        tmp_up = zeros(1,xr);
        tmp_dw = zeros(1,xr);
        tmp_tr = zeros(1,xr);
        tmp_up(trl_dw(index_y + 1,:) ~= 0) = 1;
        tmp_dw(trl_dw(index_y - 1,:) ~= 0) = 1;
        tmp_tr(trl_dw(index_y,:) ~= 0) = 1;
            if index_y == 1 || index_y == index_y_max
            index_x_b = find((tmp_up + tmp_dw + tmp_tr) == 2);
            else
            index_x_b = find((tmp_up + tmp_dw + tmp_tr) == 3);            
            end
        diff_up = zeros(1,xr);
        diff_dw = zeros(1,xr);
        diff_up(index_x_b) = trl_dw(index_y,index_x_b) -  trl_dw(index_y - 1,index_x_b);
        diff_dw(index_x_b) = trl_dw(index_y,index_x_b) -  trl_dw(index_y + 1,index_x_b);
        q_up = find(abs(diff_up) > 1.0*pi);
        q_dw = find(abs(diff_dw) > 1.0*pi);
        if isempty(q_up) || isempty(q_dw)
        else
            quality_up_tmp = (length(index_x_b) - length(q_up))/length(index_x_b);
            quality_dw_tmp = (length(index_x_b) - length(q_dw))/length(index_x_b); 
            quality_br_tmp = quality_up_tmp*quality_dw_tmp; 
        %--------------------
            bad_up = zeros(1,xr);
            bad_dw = zeros(1,xr);
            if quality_br_tmp < 0.8
            %quality_up_tmp
            %quality_dw_tmp
         %index_y
            bad_up(q_up) = 2;
            bad_dw(q_dw) = 1;         
            bad_point = bad_up + bad_dw;
%            figure(15)
%            plot(trl_dw(index_y,:),'o-r');hold on;
%            plot(bad_point,'o-b'); hold off;
%            pause
 %       index_s_dw = find(ref_dw(index_y,:) ~= 0);  
 %       %--------------------
 %       ref_dw_tmp = ref_dw(index_y,index_s_dw);
 %       trl_dw_tmp = trl_dw(index_y,index_s_dw);
 %       [index_ls_TR] = dePULM_2D_ls(trl_dw_tmp,index_s_dw); 
 %       %    
 %       [out_TR] = dePULM_2D_merging(trl_dw_tmp,index_ls_TR,ref_dw_tmp);
 %       %
 %       trl_dw(index_y,index_s_dw) = out_TR;
            end
        end
    end
%--------------
    for index_y = 2:(index_y_max - start_up + 0)
        tmp_up = zeros(1,xr);
        tmp_dw = zeros(1,xr);
        tmp_tr = zeros(1,xr);
        tmp_up(trl_up(index_y + 1,:) ~= 0) = 1;
        tmp_dw(trl_up(index_y - 1,:) ~= 0) = 1;
        tmp_tr(trl_up(index_y,:) ~= 0) = 1;
            if index_y == 1 || index_y == index_y_max
            index_x_b = find((tmp_up + tmp_dw + tmp_tr) == 2);
            else
            index_x_b = find((tmp_up + tmp_dw + tmp_tr) == 3);            
            end
        diff_up = zeros(1,xr);
        diff_dw = zeros(1,xr);
        diff_up(index_x_b) = trl_up(index_y,index_x_b) -  trl_up(index_y - 1,index_x_b);
        diff_dw(index_x_b) = trl_up(index_y,index_x_b) -  trl_up(index_y + 1,index_x_b);
        q_up = find(abs(diff_up) > 1.0*pi);
        q_dw = find(abs(diff_dw) > 1.0*pi);
        quality_up_tmp = (length(index_x_b) - length(q_up))/length(index_x_b);
        quality_dw_tmp = (length(index_x_b) - length(q_dw))/length(index_x_b); 
        quality_br_tmp = quality_up_tmp*quality_dw_tmp; 
        %--------------------
        bad_up = zeros(1,xr);
        bad_dw = zeros(1,xr);
        if quality_br_tmp < 0.8
         %index_y + start_up
 %           quality_up_tmp
 %           quality_dw_tmp         
         bad_up(q_up) = 2;
         bad_dw(q_dw) = 1;         
         bad_point = bad_up + bad_dw;
 %        figure(15)
 %        plot(trl_up(index_y,:),'o-r');hold on;
 %        plot(bad_point,'o-b'); hold off;
 %        pause
        %--------------------
 %       index_s_up = find(ref_up(index_y,:) ~= 0);
 %       %--------------------
 %       ref_up_tmp = ref_up(index_y,index_s_up);
 %       trl_up_tmp = trl_up(index_y,index_s_up);
 %       %---------------------
 %       [index_ls_TR] = dePULM_2D_ls(trl_up_tmp,index_s_up); 
 %       %    
 %       [out_TR] = dePULM_2D_merging(trl_up_tmp,index_ls_TR,ref_up_tmp);
        %
 %       trl_up(index_y,index_s_up) = out_TR;
        end
    end    
%---------------------------------------------------------------------- 
%    phase_tmp(index_y_min:start_dw,:) = trl_dw;
%    phase_tmp(start_up:index_y_max,:) = trl_up;
%-----------------
out_2D = phase_tmp;
%   phase_unwrapping was done
