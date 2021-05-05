function R2 = greedyR2(J, nVxl)
%GREEDYR2 Returns the R2* resulting in the smallest residual
%   Returns the R2* resulting in the smallest residual, assumes a single
%   minimum

nR2 = size(J,1);
R2 = zeros(nVxl,1);
for i = 1:nVxl
    r = 0;
    min = J(r+1, i);
    while r+1 < nR2 && J(r+2, i) < min
        min = J(r+2, i);
        r = r + 1;
    end
    R2(i) = r;
end

end