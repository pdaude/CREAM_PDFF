% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_2D trl_mean] = dePULM_mean_merg(phase_tmp_ref, phase_tmp,cell_signal, cell_connect,ref_mean,trl_mean,index_y_min,index_y_max)
%----------------------------------------------------------------------
    pi_2 = pi*2.0;
    pi_05 = pi*0.5;
    for index_y = index_y_min:index_y_max %1:length(index_signal)
        %   calculating the  difference of connect_mean
        diff_mean_tmp = ref_mean(index_y) - trl_mean(index_y);
        index_u = cell_signal{1,index_y};
        phase_tmp(index_y,index_u) = phase_tmp(index_y,index_u) + round(diff_mean_tmp/pi_2)*pi_2;
        index_x = cell_connect{1,index_y};
        trl_mean(index_y) = mean_liu(phase_tmp(index_y,index_x));            
    end
    %-------------------------------------------------------
    diff_mean = ref_mean - trl_mean;
    index_diff_g05pi = find(abs(ref_mean - trl_mean) > pi_05);
    if isempty(index_diff_g05pi)
    else
        for ii = 1:length(index_diff_g05pi)
        index_y = index_diff_g05pi(ii);
            phase_RF = phase_tmp_ref(index_y,:);
            phase_TR = phase_tmp(index_y,:);
            [index_ls_RF  index_RF RF_phase_tmp] = dePULM_2D_ls(phase_RF,1);
            [index_ls_TR  index_TR TR_phase_tmp] = dePULM_2D_ls(phase_TR,0); 
        %    
            ls_mag = index_TR;
            index_ls = index_ls_TR;
            line_phase = TR_phase_tmp(ls_mag);        
            re_phase = RF_phase_tmp(ls_mag);
            [out_TR] = dePULM_2D_merging(line_phase,index_ls,re_phase);
        %
        phase_TR(ls_mag) = out_TR;
        phase_tmp(index_y,:) = phase_TR;
        %------------------------
        index_x = cell_connect{1,index_y};
        trl_mean(index_y) = mean_liu(phase_tmp(index_y,index_x));   
        %------------------------
        diff_mean_tmp = ref_mean(index_y) - trl_mean(index_y);
            %if abs(diff_mean_tmp) > pi_05
            %    diff_phase_dw = phase_RF(index_x) - phase_tmp(index_y + 1,index_x);
            %    diff_phase_up = phase_RF(index_x) - phase_tmp(index_y - 1,index_x);
            %    if length(diff_phase_dw) > 0.
            %    
            %end
        end
    end
    %
    for index_y = index_y_min:index_y_max
        diff_test = phase_tmp_ref(index_y,:) - phase_tmp(index_y,:);
        index_x = cell_connect{1,index_y};
        index_diff_g10pi = find(abs(diff_test(index_x)) > pi);
        if length(index_diff_g10pi)/length(index_x) >= 1.8
        phase_RF = phase_tmp_ref(index_y,:);
        phase_TR = phase_tmp(index_y,:);
        [index_ls_RF  index_RF RF_phase_tmp] = dePULM_2D_ls(phase_RF,1);
        [index_ls_TR  index_TR TR_phase_tmp] = dePULM_2D_ls(phase_TR,0); 
        %    
        ls_mag = index_TR;
        index_ls = index_ls_TR;
        line_phase = TR_phase_tmp(ls_mag);        
        re_phase = RF_phase_tmp(ls_mag);
        [out_TR] = dePULM_2D_merging(line_phase,index_ls,re_phase);
        %
        phase_TR(ls_mag) = out_TR;
        phase_tmp(index_y,:) = phase_TR;
        %------------------------
        trl_mean(index_y) = mean_liu(phase_tmp(index_y,index_x));
        %------------------------
        end        
    end
    %--------------------------------------------------------
out_2D = phase_tmp;
%--------------------------------------------------------------------------
%   phase_unwrapping was done
