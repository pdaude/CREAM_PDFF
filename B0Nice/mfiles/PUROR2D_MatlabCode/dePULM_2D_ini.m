% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [cell_signal cell_connect cell_seg_x] = dePULM_2D_ini(mag_tmp)
%----------------------------------------------------------------------
%   load the complex imaging file
    [yr_tmp xr_tmp] = size(mag_tmp);
%----------------------------------------------------------------------
    cell_signal = cell(1,yr_tmp);
    cell_seg_x = cell(1,yr_tmp);
    for jj = 1:yr_tmp 
        index_jj = find(mag_tmp(jj,:) == 1);
        cell_signal(1,jj) = {index_jj};
        %------------------------------
        seg_x = [];
        if length(index_jj) >= 6
            kk = 0;
            end_flag = 0;
            for ii = 1:(length(index_jj) - 1)
                if index_jj(ii + 1) - index_jj(ii) == 1
                    if end_flag == 0
                    kk = kk + 1;
                    end
                    %---------
                    if kk == 3
                    seg_x = [seg_x (ii - 2)];
                    end_flag = 1;
                    kk = 0;
                    end
                    %------------
                    if ii == length(index_jj) - 1 && end_flag == 1
                    seg_x = [seg_x (ii + 1)]; 
                    end_flag = 0;
                    end
                else
                    if end_flag == 1
                    seg_x = [seg_x ii];
                    end_flag = 0;
                    end
                end
            end
        end
        cell_seg_x(1,jj) = {seg_x};
        %
    end    
    %----------------------------------------------------------------------
    cell_connect = cell(1,yr_tmp);
    for jj = 1:yr_tmp
        if jj ~= 1 && jj ~= yr_tmp
            sum_tmp_y = mag_tmp(jj+1,:) + mag_tmp(jj,:)+ mag_tmp(jj-1,:);
                bb = 1;
                mag_signal_tmp_connect = zeros(1,xr_tmp); 
            for ii = 2:(xr_tmp - 1)
                if sum_tmp_y(ii) == 3 && mag_tmp(jj,ii-1) + mag_tmp(jj,ii)+ mag_tmp(jj,ii+1) == 3
                mag_signal_tmp_connect(bb) = ii; 
                bb = bb + 1;  
                end    
            end 
        %
            if bb <= xr_tmp/4             
              for ii = 2:(xr_tmp - 1)
                if sum_tmp_y(ii) == 3
                  mag_signal_tmp_connect(bb) = ii; 
                  bb = bb + 1;   
                end              
              end
            end
            mag_signal_tmp_connect(bb:xr_tmp) = []; 
        else
             if jj == 1
             sum_tmp_y = mag_tmp(jj + 1,:) + mag_tmp(jj,:); 
             else
             sum_tmp_y = mag_tmp(jj - 1,:) + mag_tmp(jj,:); 
             end
             %
             bb = 1;            
             mag_signal_tmp_connect = zeros(1,xr_tmp);
             for ii = 2:(xr_tmp - 1)
                if sum_tmp_y(ii) == 2 && mag_tmp(jj,ii-1) + mag_tmp(jj,ii)+ mag_tmp(jj,ii+1) == 3
                mag_signal_tmp_connect(bb) = ii; 
                bb = bb + 1;  
                end    
             end 
             %
             if bb <= xr_tmp/4
              for ii = 2:(xr_tmp - 1)
                if sum_tmp_y(ii) == 2
                  mag_signal_tmp_connect(bb) = ii; 
                  bb = bb + 1;
                end              
              end
             end
            mag_signal_tmp_connect(bb:xr_tmp) = []; 
        end
        %
        if bb ~= 1
        cell_connect(1,jj) = {mag_signal_tmp_connect};
        end
    end
    %
%--------------------------------------------------------------------

%--------------------------------------------------------   
%out_2D = unwrapped_phase;
%--------------------------------------------------------------------------
%   Initializing phase_2D was done
