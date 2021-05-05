% This program is to stack the multiple slices unwrapped in 2D into a 3D
% volume by unwrapping the mean values of the central line and the shift 
%the phase values of the volume
function [line_unwrapping cell_mean] = dePULM_2D_stacking(line_unwrapping, line_mag, cell_mean,slice_start, slice_end,base_line)
%
pi_2 = 2.0*pi;
total_slice = slice_end - slice_start + 1;
%--------------
D3_mean = zeros(1,total_slice);
for index_slice = slice_start:slice_end
    D3_mean_test = cell_mean{1,index_slice}(base_line);
    D3_mean_test = D3_mean_test + cell_mean{1,index_slice}(base_line - 1);
    D3_mean_test = D3_mean_test + cell_mean{1,index_slice}(base_line + 1);    
    D3_mean(index_slice) = D3_mean_test/3;%cell_mean{1,index_slice}(base_line);
end
%
unwrap_D3_mean = unwrap(D3_mean);
central_slice = round((slice_start + slice_end)/2);
unwrap_D3_mean = unwrap_D3_mean - pi_2*round(unwrap_D3_mean(central_slice)/(pi_2));
%--------------
for index_slice = 1:total_slice
    diff_test = unwrap_D3_mean(index_slice) - D3_mean(index_slice);
    pi_w = pi;
    if abs(diff_test) > pi_w
    line_unwrapping{index_slice}(:,:) = line_unwrapping{index_slice}(:,:).*line_mag{index_slice}(:,:)...
                                        + pi*round(diff_test/(pi)).*line_mag{index_slice}(:,:);
    cell_mean_tmp = cell_mean{1,index_slice} + pi*round(diff_test/(pi));
    index_cell_mean_tmp = find(cell_mean_tmp == pi*round(diff_test/(pi)));
    cell_mean_tmp(index_cell_mean_tmp) = 0;
    cell_mean{1,index_slice} = cell_mean_tmp;
    end
end 
    %--------------------------------------

    