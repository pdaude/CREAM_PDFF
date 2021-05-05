% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [index_ls] = dePULM_2D_ls(phase_ls,index_u,phi_good)
%----------------------------------------------------------------------
            phase_1D = phase_ls;
            %phi_good = pi*1.0;
            w_good = 2;
            %
            length_1D = length(phase_1D);
            %
if length_1D >= 6
            %--------------------------
            index_ls_tmp = zeros(1,length_1D);
            kk = 1;
            for ii = 1:(length(phase_1D)-1)
                diff_nb = phase_1D(ii+1) - phase_1D(ii);
                if abs(diff_nb) <= (phi_good*1.0) && (index_u(ii+1) - index_u(ii) == 1)
                   index_ls_tmp(kk) = ii;
                   kk = kk + 1;
                end
            end
            index_ls_tmp(kk:length_1D) = [];
            %
        if isempty(index_ls_tmp)
        index_ls = [1 length(phase_1D)];
        else
            index_ls_tmp_tmp = zeros(1, length_1D);
            index_ls_tmp_tmp(1) = index_ls_tmp(1);
            kk = 2;
            for ii = 2:length(index_ls_tmp)
                if (index_ls_tmp(ii) - index_ls_tmp(ii - 1) > 1)
                    index_ls_tmp_tmp(kk) = (index_ls_tmp(ii-1) + 1);
                    index_ls_tmp_tmp(kk+1) = (index_ls_tmp(ii) + 0);
                    kk = kk + 2;
                end    
            end
            index_ls_tmp_tmp(kk) = index_ls_tmp(length(index_ls_tmp));
            index_ls_tmp_tmp((kk+1):length_1D) = [];
            %
            index_ls = zeros(1, length_1D);
            kk = 1;
            %
            if length(index_ls_tmp_tmp) > w_good
                for ii = 1:2:length(index_ls_tmp_tmp)
                    if index_ls_tmp_tmp(ii + 1) - index_ls_tmp_tmp(ii) > 2
                       index_ls(kk:(kk+1))  = index_ls_tmp_tmp(ii:(ii+1)); 
                       kk = kk + 2; 
                    end    
                end
                index_ls(kk:length_1D) = [];
            else
                index_ls = index_ls_tmp_tmp;
            end
        end
%-------------------------------------------------        
else
    index_ls = [1 length(phase_1D)];
end
    %
if isempty(index_ls)
    index_ls = [1 length(phase_1D)];
end
%
%------------------------ 


