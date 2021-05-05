function J = getR2Residuals(Y, dB0, C, nB0, nR2, nVxl)
%GETR2RESIDUALS Gets R2* residuals
%   Gets R2* residuals

J = zeros(nR2, nVxl);
for b = 1:nB0
    for r = 1:nR2
        y = Y(:, dB0 == b);
        J(r, dB0 == b) = sum(abs((squeeze(C(r,b,:,:)) * y)).^2);
    end
end

end