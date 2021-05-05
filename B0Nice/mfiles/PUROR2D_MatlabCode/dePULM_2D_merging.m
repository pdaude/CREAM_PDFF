function [line_tmp_phase] = dePULM_2D_merging(line_phase,index_ls,re_phase) 
%
    pi_2 = pi*2.0;
    phi_g = 1.0*pi;
    line_tmp_phase = line_phase;
if length(index_ls) > 2
%-----------
    diff_tmp = re_phase - line_phase;
%------------------------------------------------------------------------
    seg_mean_shift = zeros(1,length(index_ls));
    %
    for ii = 1:2:(length(index_ls) - 1)
    kk = 0;
    sum_diff = 0;
        for jj = index_ls(ii):index_ls(ii+1)
            %if re_phase(jj) ~= 0
                sum_diff = sum_diff + diff_tmp(jj);
                kk = kk + 1;
            %end
        end    
    % 
        if kk > 0
            ave_diff = sum_diff/kk;
            seg_mean_shift(ii) = ave_diff;
            seg_mean_shift(ii+1) = ave_diff;           
        end
    end
%-------------------------------------------------------
  if isempty(abs(seg_mean_shift) > phi_g) 
  else
    %
    for ii = 1:2:(length(index_ls) - 1)
       ave_diff = seg_mean_shift(ii);
       ls_s = index_ls(ii);
       ls_e = index_ls(ii+1);
       line_tmp_phase(ls_s:ls_e) = line_phase(ls_s:ls_e) + pi_2*round(ave_diff/(pi_2));  
    end
    %
    for ii = 2:2:(length(index_ls) - 2)
       in_s = index_ls(ii) + 1;
       in_e = index_ls(ii+1) -1;
       if in_e - in_s >= 0
           kk = 0; 
           sum_diff = 0;
            for jj = (index_ls(ii)-1):(index_ls(ii+1)+1);
                %if re_phase(jj) ~= 0
                 sum_diff = sum_diff + diff_tmp(jj);
                 kk = kk + 1;
                %end
            end
            if kk == 0
            ave_diff = (seg_mean_shift(ii) + seg_mean_shift(ii+1))/2;
            else
            ave_diff =  sum_diff/kk;%(seg_mean_shift(ii) + seg_mean_shift(ii+1))/2;
            end
            line_tmp_phase(in_s:in_e) = line_phase(in_s:in_e) + pi_2*round(ave_diff/(pi_2));
       end
    end
    %    
    %
  end%
  %
end
%--------------------------------------------------------------------------