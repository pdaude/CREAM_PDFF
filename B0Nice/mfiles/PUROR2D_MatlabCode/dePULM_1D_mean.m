% this program is for volunteer results acquired using 3D SPGR in GE 3T
% scanner
%
%
function [mean_connect_y mean_connect_x] = dePULM_1D_mean(phase_tmp,cell_connect_x,cell_connect_y)
%
[yr_tmp xr_tmp] = size(phase_tmp);
%----------------------------------------------------------------------
    mean_connect_y = zeros(1,yr_tmp);
    for index_y = 1:yr_tmp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   calculating the connect_mean
        index_c = cell_connect_x{1,index_y};
        if isempty(index_c)
        else
        phase_x = phase_tmp(index_y,:);
        mean_connect_y(index_y) = mean_liu(phase_x(index_c)); 
        end
    end
    %--------------------------------------------------------
    mean_connect_x = zeros(1,yr_tmp);
    for index_x = 1:xr_tmp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   calculating the connect_mean
        index_c = cell_connect_y{1,index_x};
        if isempty(index_c)
        else
        phase_y = phase_tmp(:,index_x);
        mean_connect_x(index_x) = mean_liu(phase_y(index_c));
        end
    end    