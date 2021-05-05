% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [cell_signal cell_connect cell_seg_y] = dePULM_2D_ini_y(mag_tmp)
%----------------------------------------------------------------------
%   load the complex imaging file
    [yr_tmp xr_tmp] = size(mag_tmp);
%----------------------------------------------------------------------
    cell_signal = cell(1,xr_tmp);
    cell_seg_y = cell(1,xr_tmp);
    for index_x = 1:xr_tmp 
        index_jj = find(mag_tmp(:,index_x) == 1);
        cell_signal(1,index_x) = {index_jj};
        %------------------------------
        seg_y = [];
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
                    seg_y = [seg_y (ii - 2)];
                    end_flag = 1;
                    kk = 0;
                    end
                    %------------
                    if ii == length(index_jj) - 1 && end_flag == 1
                    seg_y = [seg_y (ii + 1)]; 
                    end_flag = 0;
                    end
                else
                    if end_flag == 1
                    seg_y = [seg_y ii];
                    end_flag = 0;
                    end
                end
            end
        end
        cell_seg_y(1,index_x) = {seg_y};
        %
    end    
    %----------------------------------------------------------------------
    cell_connect = cell(1,xr_tmp);
    for index_x = 1:xr_tmp
        if index_x ~= 1 && index_x ~= xr_tmp
        sum_tmp_x = mag_tmp(:,index_x + 1) + mag_tmp(:,index_x)+ mag_tmp(:,index_x - 1);
        mag_signal_tmp_connect = zeros(1,yr_tmp);
        bb = 1;
        for ii = 2:(yr_tmp - 1)
            if sum_tmp_x(ii) == 3 && mag_tmp(ii-1,index_x) + mag_tmp(ii,index_x)+ mag_tmp(ii+1,index_x) == 3
                  mag_signal_tmp_connect(bb) = ii; 
                  bb = bb + 1;  
            end    
        end
        %      
        if bb <= yr_tmp/4
              for ii = 2:(yr_tmp - 1)
                if sum_tmp_x(ii) == 3
                  mag_signal_tmp_connect(bb) = ii; 
                  bb = bb + 1;
               %end    
                end              
              end
        end
        mag_signal_tmp_connect(bb:yr_tmp) = []; 
        else
            if index_x == 1
            sum_tmp_x = mag_tmp(:,index_x + 1) + mag_tmp(:,index_x);           
            else
            sum_tmp_x = mag_tmp(:,index_x)+ mag_tmp(:,index_x - 1);    
            end
            mag_signal_tmp_connect = zeros(1,yr_tmp);
            bb = 1;
            for ii = 2:(yr_tmp - 1)
                if sum_tmp_x(ii) == 2 && mag_tmp(ii-1,index_x) + mag_tmp(ii,index_x)+ mag_tmp(ii+1,index_x) == 3
                  mag_signal_tmp_connect(bb) = ii; 
                  bb = bb + 1;  
                end    
            end
        %      
            if bb <= yr_tmp/4
              for ii = 2:(yr_tmp - 1)
                if sum_tmp_x(ii) == 2
                  mag_signal_tmp_connect(bb) = ii; 
                  bb = bb + 1;
               %end    
                end              
              end
            end
            mag_signal_tmp_connect(bb:yr_tmp) = [];                    
        end
        %
        if bb ~= 1
        cell_connect(1,index_x) = {mag_signal_tmp_connect};
        end
    end
    %--------------------------------------------------------   
%out_2D = unwrapped_phase;
%--------------------------------------------------------------------------
%   Initializing phase_2D was done
