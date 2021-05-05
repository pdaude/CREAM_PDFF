% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [out_2D] = dePULM_2D_mean(phase_original,cell_signal, cell_connect, cell_seg_x)
%----------------------------------------------------------------------
    pi_2 = pi*2.0;
    [yr_tmp xr_tmp] = size(phase_original);
    cen_line = round(yr_tmp/2);
    %---------------------------------------------------------
        phase_tmp = zeros(yr_tmp,xr_tmp);
    for index_y = 1:yr_tmp
        index_u = cell_signal{1,index_y};
        len_u = length(index_u);
        if len_u > xr_tmp/4
        phase_u = phase_original(index_y,index_u);
        g_seg = cell_seg_x{1,index_y};
        phi_good = 0.5*pi;
        [index_ls] = dePULM_1D_ls(phase_u,phi_good,g_seg);
        [phase_tmp_1] = dePULM_1D(phase_u,index_u,index_ls);
        phase_tmp(index_y,index_u) = phase_tmp_1(1:len_u);
        else
        phase_tmp(index_y,index_u) = phase_original(index_y,index_u);            
        end
    end
    %-------------------------------------------------------
    mean_connect = zeros(1,yr_tmp);
    for index_y = 2:yr_tmp %1:length(index_signal)
        phase_x = phase_tmp(index_y,:);
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   calculating the connect_mean
        index_c = cell_connect{1,index_y};
        index_s = cell_signal{1,index_y};
        if isempty(index_s)
        mean_connect(index_y) = mean_connect(index_y - 1);        
        else
        %
            if isempty(index_c) || length(index_c) <= 3
            mean_connect(index_y) = mean_liu(phase_x(index_s));        
            else
            mean_connect(index_y) = mean_liu(phase_x(index_c));
            end
        %
        phase_tmp(index_y,index_s) = phase_tmp(index_y,index_s) - round(mean_connect(index_y)/pi_2)*pi_2;
        %
        phase_x = phase_tmp(index_y,:);
            if isempty(index_c) || length(index_c) <= 3
            mean_connect(index_y) = mean_liu(phase_x(index_s));        
            else
            mean_connect(index_y) = mean_liu(phase_x(index_c));
            end
        %
        end
    end
    %--------------------------------------------------------
    index_y = 1;
        phase_x = phase_tmp(index_y,:);
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   calculating the connect_mean
        index_c = cell_connect{1,index_y};
        index_s = cell_signal{1,index_y};
        if isempty(index_s)
        mean_connect(index_y) = mean_connect(index_y + 1);        
        else
        %
            if isempty(index_c) || length(index_c) <= 3
            mean_connect(index_y) = mean_liu(phase_x(index_s));        
            else
            mean_connect(index_y) = mean_liu(phase_x(index_c));
            end
        %
        phase_tmp(index_y,index_s) = phase_tmp(index_y,index_s) - round(mean_connect(index_y)/pi_2)*pi_2;
        %
        phase_x = phase_tmp(index_y,:);
            if isempty(index_c) || length(index_c) <= 3
            mean_connect(index_y) = mean_liu(phase_x(index_s));        
            else
            mean_connect(index_y) = mean_liu(phase_x(index_c));
            end
        %
        end    
%--------------------------------------------------------------------
%   Start line shift by unwrapping the means    
%--------------------------------------------------------------------
%   unwrap the mean values before global shfit the phase data
    mean_u = mean_connect(1:yr_tmp);
    mean_unwrap = zeros(1,yr_tmp);
    mean_unwrap(1:yr_tmp) = unwrap(mean_u);
%----------------------------------------------------------
    mean_unwrap(1:yr_tmp) = mean_unwrap(1:yr_tmp)...
                                         - round(mean_unwrap(cen_line)/pi_2)*pi_2;
%--------------------------------------------------------------------------
%   shift the phase data
    for index_y = 1:yr_tmp
        index_s = cell_signal{1,index_y};        
        diff_test = mean_unwrap(index_y) - mean_connect(index_y);
        if abs(diff_test) > pi
        phase_tmp(index_y,index_s) = phase_tmp(index_y,index_s)...
                               + pi_2*round(diff_test/(pi_2));                                      
        end
        %
    end
%
out_2D = phase_tmp;
%--------------------------------------------------------------------------
%   phase_unwrapping was done
