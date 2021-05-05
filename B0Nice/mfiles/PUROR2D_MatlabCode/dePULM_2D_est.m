% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_2D mean_tmp] = dePULM_2D_est(unwrapped_phase, mag_tmp,cell_connect, mean_tmp,index_y_min,index_y_max)
%----------------------------------------------------------------------
%   load the imaging file
    pi_2 = pi*2;
    phase_range = 4.0*pi; %[-phase_range phase_range]
    [yr xr] = size(unwrapped_phase);
    jump_w = pi;
%--------------------------------------------------------------
   for jj_loop = 1:2 
        %
        if mod(jj_loop,2) == 1
        [q_br_y q_up_y q_dw_y gd_y] = dePULM_2D_quality(unwrapped_phase,index_y_min,index_y_max);        
        protect_start = gd_y(1);
        protect_end = gd_y(2);
        protect_start_tmp = gd_y(3);
        protect_end_tmp = gd_y(4);
        end
        %
        if mod(jj_loop,2) == 1
        nxt_drn = -1;
        ref_drn = +1;
        jj = protect_start_tmp;
        else
        nxt_drn = +1;
        ref_drn = -1;
        jj = protect_end_tmp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        do_2D = 1;
        while do_2D == 1
            phase_RF = unwrapped_phase(jj + ref_drn,:);
            phase_TR = unwrapped_phase(jj,:);
            phase_NT = unwrapped_phase(jj + nxt_drn,:);
            flag_display = 0;
            %           
            if mod(jj_loop,2) == 1
            trl_ref = q_dw_y(jj);
            else
            trl_ref = q_up_y(jj);               
            end
            q_trl = q_br_y(jj);
            %
            if jj_loop <= 2
                th_trl = 0.8;
            else
                th_trl = 0.8;
            end
            %
        if  trl_ref <= th_trl
            if flag_display == 1
            display ('need 2D');
            display ('step 0');
            figure(9)
            plot(phase_RF,'o-k'); hold on;
            plot(phase_TR,'o-r'); hold on;
            plot(phase_NT,'o-g'); hold off;
            pause
            end
            %--------------------------------------------------------------
            if jj_loop <= 2
            [index_ls_TR  index_TR TR_phase_tmp] = dePULM_2D_ls(phase_TR,0); 
            %
            ls_mag = index_TR;
            index_ls = index_ls_TR;
            line_phase = TR_phase_tmp(ls_mag);        
            re_phase = phase_RF(ls_mag);
            [out_TR] = dePULM_2D_merging(line_phase,index_ls,re_phase);
            %
            tmp_line_phase = zeros(1,xr);
            tmp_line_phase(ls_mag) = out_TR;
   %---------------------------------------------------------------
            %       Method 2
         else
            index_jj = cell_connect{1,jj};
            if isempty(index_jj)
            else
          %  jump_w = pi*1.0;
            gg_w = 24;
            tmp_tmp = zeros(1,xr);
            tmp_tmp(index_jj) = 1;
            mag_jj = mag_tmp(jj,:);
%            phase_TR = tmp_line_phase; % using both two methods
            diff_test = (phase_RF - phase_TR).*tmp_tmp;
            jump_tmp  = find(abs(diff_test) > jump_w);
            %--------------------------------------------------------------
                if isempty(jump_tmp)
                tmp_line_phase = phase_TR;
                else
                [phase_jj_out] = dePU_line_merging(phase_TR,mag_jj,index_jj,diff_test,jump_tmp,gg_w,xr);
                tmp_line_phase = phase_jj_out;
                end
            end    
          end
            %recaculated the quality;
            index_jj = cell_connect{1,jj};
            if isempty(index_jj)
            else
            ref_phase = unwrapped_phase(jj+ref_drn,index_jj);
            trl_phase = tmp_line_phase(index_jj);
            nxt_phase = unwrapped_phase(jj + nxt_drn,index_jj);
            diff_nxt = trl_phase - nxt_phase;
            diff_ref = trl_phase - ref_phase;
            trl_nxt = length(find(abs(diff_nxt) < 1.0*pi))/length(index_jj);
            trl_ref = length(find(abs(diff_ref) < 1.0*pi))/length(index_jj);
                if mod(jj_loop,2) == 1
                q_up_y(jj) = trl_nxt;
                q_dw_y(jj-1) = trl_nxt;
                q_dw_y(jj) = trl_ref;
                q_br_y(jj) = trl_nxt*trl_ref;
                q_br_y(jj - 1) = q_dw_y(jj-1)*q_up_y(jj-1);
                else
                q_up_y(jj) = trl_ref;
                q_dw_y(jj) = trl_nxt;
                q_up_y(jj+1) = trl_nxt;
                q_br_y(jj) = trl_nxt*trl_ref;
                q_br_y(jj + 1) = q_dw_y(jj + 1)*q_up_y(jj + 1);            
                end
            %--------------------------------------------------------------
                if flag_display == 1
                display ('method 1');
                figure(10)
                plot(ref_phase,'o-k'); hold on;
                plot(trl_phase,'o-r'); hold on;
                plot(nxt_phase,'o-g'); hold off;           
                pause
                end
            %------------------------------------------------------------- 
            unwrapped_phase(jj,:) = tmp_line_phase;
            end
        end
            %
            index_jj = cell_connect{1,jj};
            jj_tmp = unwrapped_phase(jj,index_jj);
            mean_tmp(jj) = mean_liu(jj_tmp);  
            %
            if mod(jj_loop,2) ==1
                if jj > (index_y_min + 2)
                mean_diff = mean_tmp(jj + ref_drn) - mean_tmp(jj);
                    if abs(mean_diff) > jump_w*10000
                    kk = jj;
                    unwrapped_phase(index_y_min:kk,:) = unwrapped_phase(index_y_min:kk,:)...
                               + pi_2*round(mean_diff/(pi_2)).*mag_tmp(index_y_min:kk,:);
                    mean_tmp(index_y_min:kk) = mean_tmp(index_y_min:kk)...
                               + pi_2*round(mean_diff/(pi_2));
                    end
                end
                            %
               if jj > (index_y_min + 2)
               mean_diff = mean_tmp(jj) - mean_tmp(jj +nxt_drn);
                    if abs(mean_diff) > jump_w
                    kk = jj + nxt_drn;
                    unwrapped_phase(index_y_min:kk,:) = unwrapped_phase(index_y_min:kk,:)...
                               + pi_2*round(mean_diff/(pi_2)).*mag_tmp(index_y_min:kk,:);
                    mean_tmp(index_y_min:kk) = mean_tmp(index_y_min:kk)...
                               + pi_2*round(mean_diff/(pi_2));
                    end
               end
               %
            else
                %
                if jj < (index_y_max - 2)
                mean_diff = mean_tmp(jj + ref_drn) - mean_tmp(jj);
                    if abs(mean_diff) > jump_w*100000
                    kk = jj;
                unwrapped_phase(kk:index_y_max,:) = unwrapped_phase(kk:index_y_max,:)...
                               + pi_2*round(mean_diff/(pi_2)).*mag_tmp(kk:index_y_max,:);
                mean_tmp(kk:index_y_max) = mean_tmp(kk:index_y_max)...
                               + pi_2*round(mean_diff/(pi_2));
                    end
                end
                %
                if jj < (index_y_max - 2)
                mean_diff = mean_tmp(jj) - mean_tmp(jj + nxt_drn);
                    if abs(mean_diff) > jump_w*1.3
                    kk = jj + nxt_drn;
                    unwrapped_phase(kk:index_y_max,:) = unwrapped_phase(kk:index_y_max,:)...
                               + pi_2*round(mean_diff/(pi_2)).*mag_tmp(kk:index_y_max,:);
                    mean_tmp(kk:index_y_max) = mean_tmp(kk:index_y_max)...
                               + pi_2*round(mean_diff/(pi_2));
                    end
                end
            end
            %--------------------------------------------------------------
            jj = jj + nxt_drn;
            if jj < index_y_min
            do_2D = 0;
            end
            if jj > index_y_max
            do_2D = 0;
            end
            %  
        end     
   end
%    figure(41)
%    plot(mean_tmp/pi,'o-r');  hold off
%    plot(diff(mean_tmp)/pi,'o-k');  hold off;
    %-------------------------------------------------
%--------------------------------------------------------------------
%   For the final step, Start line shift by unwrapping the means    
%-------------------------------------------------------------------- 
out_2D = unwrapped_phase;
%--------------------------------------------------------------------------
%   phase_unwrapping was done
