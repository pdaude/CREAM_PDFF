% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [index_ls] = dePULM_1D_ls(phase_1D,phi_good,g_seg)
%----------------------------------------------------------------------
            %phase_1D = phase_ls;
            w_good = 3;
            %
            length_1D = length(phase_1D);
            %
if isempty(g_seg)
    if isempty(phase_1D)
    index_ls = [];
    else
    index_ls = [1 length(phase_1D)];        
    end
else
            %--------------------------
            index_ls = [];
            for jj = 1:2:(length(g_seg)-1)
                    kk = 0;
                    end_flag = 0;
                for ii = g_seg(jj):(g_seg(jj+1) - 1)
                    diff_nb = phase_1D(ii+1) - phase_1D(ii);
                    if abs(diff_nb) <= (phi_good*1.0) && end_flag == 0
                    kk = kk + 1;
                    end
                    %
                    if kk == w_good
                    index_ls = [index_ls (ii - 2)];
                    kk = 0;
                    end_flag = 1;
                    end
                    %
                    if end_flag == 1 && abs(diff_nb) > phi_good
                    index_ls = [index_ls ii];
                    end_flag = 0;
                    end
                    %
                    if end_flag == 1 && ii == (g_seg(jj+1) - 1)
                    index_ls = [index_ls (ii + 1)];                   
                    end
                    
                end
            end
            %
            if isempty(index_ls)
            index_ls = [1 length(phase_1D)];
            end
%-------------------------------------------------        
end
%
%------------------------ 


