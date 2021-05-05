function [A, B] = findTwoSmallestMinima(J)
%FINDTWOSMALLESTMINIMA Finds the two smallest minima
%   Returns the two indexes corresponding to the two off-resonance values
%   with the locally (as seen over the field map space) smallest residuals 
%   for each voxel

nVxl = size(J,2);
A = zeros(nVxl, 1);
B = zeros(nVxl, 1);
for i = 1:nVxl
    f = J(:,i);
    I = findMinima(f);
    [~,II] = sort(f(I));
    minima = I(II);
    if length(minima) >= 2
        A(i) = minima(1);
        B(i) = minima(2);
    elseif length(minima) == 1
        A(i) = minima;
        B(i) = minima;
    else
        A(i) = 1;
        B(i) = 1;
    end
end

end