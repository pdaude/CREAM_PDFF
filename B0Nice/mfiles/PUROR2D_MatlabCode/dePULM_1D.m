function [phase_1D] = dePULM_1D(phase_1D,index_u,index_ls)
%
            pi_2 = 2*pi;
            phase_1D_tmp = phase_1D;
            %--------------------------
            jj = 10000;
     if length(index_ls) > 2
           for ii = 2:2:(length(index_ls) - 2)
                   kk = index_ls(ii + 1) - index_ls(ii);     
                   %55555555555555
               if kk <= 6
                     test = phase_1D(index_ls(ii + 1)) + phase_1D(index_ls(ii + 1) + 1) + phase_1D(index_ls(ii + 1) + 2)...
                          - phase_1D(index_ls(ii)) - phase_1D(index_ls(ii) - 1) - phase_1D(index_ls(ii) - 2);
                     test = test/3; 
                   if abs(test) > pi
                   jj = index_ls(ii) + 2;
                   phase_1D(jj:length(phase_1D)) = phase_1D(jj:length(phase_1D)) - round(test/pi_2)*pi_2;
                   end                    
               end                   
                   %55555555555555
               if jj ~= 10000
               inter_diff_l = phase_1D(jj) - phase_1D(jj - 1);
                   if abs(inter_diff_l) >= pi_2
                   phase_1D(jj) = phase_1D(jj) - pi_2*fix(inter_diff_l/pi_2);
                   end
               end
            end  
            %-----
            tmp_1D = index_ls(1);
            if tmp_1D ~= 1;
                index_ls = [1 (tmp_1D - 1) index_ls];
                if tmp_1D == 2
                diff_1D = phase_1D(tmp_1D) - phase_1D_tmp(1);
                phase_1D(1) = phase_1D_tmp(1) + pi_2*round(diff_1D/pi_2);        
                else
                tmp_smooth = phase_1D_tmp(1:(tmp_1D - 1));
                diff_1D = phase_1D(tmp_1D) - mean_liu(tmp_smooth);
                phase_1D(1:(tmp_1D - 1)) = tmp_smooth + pi_2*round(diff_1D/pi_2);
                end
            end
%   
           tmp_1D = index_ls(length(index_ls));
           if tmp_1D ~= length(phase_1D)
               if tmp_1D == length(phase_1D) - 1
               diff_1D = phase_1D(tmp_1D) - phase_1D_tmp(tmp_1D + 1);
               phase_1D(tmp_1D + 1) = phase_1D_tmp(tmp_1D + 1) + pi_2*round(diff_1D/pi_2);         
               else
               tmp_smooth = phase_1D_tmp((tmp_1D + 1):length(phase_1D));               
               diff_1D = phase_1D(tmp_1D) - mean_liu(tmp_smooth);
               phase_1D((tmp_1D+1):length(phase_1D)) = phase_1D_tmp((tmp_1D+1):length(phase_1D)) + pi_2*round(diff_1D/pi_2);
               end
           end
     end
%-------------------------------------------------        
%------------------------

