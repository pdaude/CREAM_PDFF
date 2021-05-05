function [mean_out] = mean_liu(data)
% this program just for the calculation of the mean values
num_points = length(data);
sum_tmp = 0;
for ii = 1:num_points
    if isnan(data(ii))
    data(ii) = 0;
    end
    sum_tmp = sum_tmp + data(ii);
end
mean_out = sum_tmp/num_points;
if isnan(mean_out)
    mean_out = 0;
end    