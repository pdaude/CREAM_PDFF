% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_2D] = dePULM_refY_trlX_final(phase_tmp_ref, phase_tmp,seg_phi,cell_seg_x,sig_x)
%----------------------------------------------------------------------
    [yr_tmp xr_tmp] = size(phase_tmp);
%----------------------------------------------------------------------
    for index_y = 1:yr_tmp  
        %--------------------
        g_seg = cell_seg_x{1,index_y};
        if isempty(g_seg)
        %-------------------
        %if isempty(phase_tmp_ref(index_y,:) ~= 0) 
        else
        index_s_dw = sig_x{1,index_y};
        %index_s_dw = find(phase_tmp_ref(index_y,:) ~= 0);
        ref_tmp = phase_tmp_ref(index_y,index_s_dw);
        trl_tmp = phase_tmp(index_y,index_s_dw);
        diff_tmp = ref_tmp - trl_tmp;
        diff_index = diff_tmp(abs(diff_tmp) > pi);
            if isempty(diff_index) || length(diff_index) < 3
            else
            [index_ls_TR] = dePULM_1D_ls(trl_tmp,seg_phi,g_seg);
            %[index_ls_TR] = dePULM_2D_ls(trl_tmp,index_s_dw,seg_phi);
            %    
            [out_TR] = dePULM_2D_merging(trl_tmp,index_ls_TR,ref_tmp);
            %
            phase_tmp(index_y,index_s_dw) = out_TR;
            end
        end
    end
%-----------------
out_2D = phase_tmp;
%   phase_unwrapping was done
