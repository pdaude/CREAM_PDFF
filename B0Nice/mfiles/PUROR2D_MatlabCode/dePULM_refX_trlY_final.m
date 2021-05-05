% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_2D] = dePULM_refX_trlY_final(phase_tmp_ref, phase_tmp,seg_phi,seg_y,sig_y)
%----------------------------------------------------------------------
    [yr_tmp xr_tmp] = size(phase_tmp);
%----------------------------------------------------------------------   
    for index_x = 1:xr_tmp
        g_seg = seg_y{1,index_x};
        if isempty(g_seg) 
        else
        %--------------------
        index_s_dw = sig_y{1,index_x};
        %index_s_dw = find(phase_tmp_ref(:,index_x) ~= 0);
        ref_tmp = phase_tmp_ref(index_s_dw,index_x);
        trl_tmp = phase_tmp(index_s_dw,index_x);
        %
        [index_ls_TR] = dePULM_1D_ls(trl_tmp,seg_phi,g_seg);
        %[index_ls_TR] = dePULM_2D_ls(trl_tmp,index_s_dw,seg_phi); 
        %    
        [out_TR] = dePULM_2D_merging(trl_tmp,index_ls_TR,ref_tmp);
        %
        phase_tmp(index_s_dw,index_x) = out_TR;
        end
    end    
%---------------------------------------------------------------------- 
out_2D = phase_tmp;
%   phase_unwrapping was done
