% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [mag_tmp MASK_phase] = dePULM_2D_mask(filter_mag, filter_th,phase_original)
%----------------------------------------------------------------------
%   load the complex imaging file
    [yr_tmp xr_tmp] = size(filter_mag);
%----------------------------------------------------------------------
    mag_tmp = zeros(1,(yr_tmp*xr_tmp));
    for ii = 1:(yr_tmp*xr_tmp)
        %if (mag_original(ii) >= thresh_mag) && (filter_mag(ii) >= filter_th)
        if (filter_mag(ii) >= filter_th)
        mag_tmp(ii) = 1; 
        end       
    end
    %
    mag_tmp = reshape(mag_tmp,yr_tmp,xr_tmp);
    %----------------------------------------------------
    mag_tmp_x = zeros(yr_tmp,xr_tmp);
    for jj = 1:yr_tmp
        if jj ~= 1 && jj ~= yr_tmp
        sum_tmp_y = mag_tmp(jj+1,:) + mag_tmp(jj,:)+ mag_tmp(jj-1,:);       
            for ii = 2:(xr_tmp - 1)
                if sum_tmp_y(ii) >= 2 && mag_tmp(jj,ii-1) + mag_tmp(jj,ii)+ mag_tmp(jj,ii+1) >= 2
               mag_tmp_x(jj,ii) = 1;
                end       
            end
        else
            sum_tmp_y = mag_tmp(jj,:);                  
            for ii = 2:(xr_tmp - 1)
                if sum_tmp_y(ii) == 1 && mag_tmp(jj,ii-1) + mag_tmp(jj,ii)+ mag_tmp(jj,ii+1) >= 2
                mag_tmp_x(jj,ii) = 1;
                end       
            end
        end
        %
        if mag_tmp_x(jj,2) == 1
           mag_tmp_x(jj,1) = 1;
        end
        %
        if mag_tmp_x(jj,(xr_tmp-1)) == 1
           mag_tmp_x(jj,xr_tmp) = 1;
        end
        %
    end    
    %----------------------------------------------------------------
    mag_tmp = rot90(mag_tmp,1);
   [yr_tmp xr_tmp] = size(mag_tmp);
    %----------------------------------------------------
    mag_tmp_y = zeros(yr_tmp,xr_tmp);
    for jj = 1:yr_tmp
        if jj ~= 1 && jj ~= yr_tmp
        sum_tmp_y = mag_tmp(jj+1,:) + mag_tmp(jj,:)+ mag_tmp(jj-1,:);       
            for ii = 2:(xr_tmp - 1)
                if sum_tmp_y(ii) >= 2 && mag_tmp(jj,ii-1) + mag_tmp(jj,ii)+ mag_tmp(jj,ii+1) >= 2
               mag_tmp_y(jj,ii) = 1;
                end       
            end
        else
        sum_tmp_y = mag_tmp(jj,:);
            for ii = 2:(xr_tmp - 1)
                if sum_tmp_y(ii) == 1 && mag_tmp(jj,ii-1) + mag_tmp(jj,ii)+ mag_tmp(jj,ii+1) >= 2
               mag_tmp_y(jj,ii) = 1;
                end       
            end        
        end
        %
        if mag_tmp_y(jj,(xr_tmp-1)) == 1
           mag_tmp_y(jj,xr_tmp) = 1;
        end        
        %
         if mag_tmp_y(jj,(xr_tmp-1)) == 1
           mag_tmp_y(jj,xr_tmp) = 1;
        end       
        %
    end
    mag_tmp_y = rot90(mag_tmp_y,3);
    [yr_tmp xr_tmp] = size(mag_tmp_y);
    %------------------------------------------------
    mag_tmp = zeros(yr_tmp,xr_tmp);
    MASK_phase = zeros(yr_tmp,xr_tmp);
    for jj = 1:yr_tmp    
        mag_tmp(jj,(mag_tmp_x(jj,:) + mag_tmp_y(jj,:) >= 1)) = 1;
       % index_jj = find(mag_tmp(jj,:) == 1);
        MASK_phase(jj,mag_tmp(jj,:) == 1) = phase_original(jj,mag_tmp(jj,:) == 1);
    end    
%--------------------------------------------------------------------------
%   Initializing phase_2D was done
