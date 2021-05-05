function J = getB0Residuals(Y, B, Bh, nB0, nVxl, R2cand, N, M, t1, dt, gyro, B0, alpha, CS, P)
%GETB0RESIDUALS Gets B0 Residuals
%   Gets B0 Residuals

RA = zeros(1,N,M);
RAp = zeros(1,M,N);
R2 = R2cand;
RA(1,:,:) = modelMatrix(N, t1, dt, M, R2, gyro, B0, alpha, CS, P);
RAp(1,:,:) = pinv(squeeze(RA(1,:,:)));

C = zeros(1,nB0,N,N);
% Null space projection matrix
proj = eye(N) - squeeze(RA(1,:,:)) * squeeze(RAp(1,:,:));
for b = 0:nB0-1
    C(1,b+1,:,:) = squeeze(B(b+1,:,:)) * proj * squeeze(Bh(b+1,:,:));
end

J = zeros(nB0, nVxl);
for b = 0:nB0-1
    J(b+1, :) = sum(abs(squeeze(C(1,b+1,:,:)) * Y).^2);
end

end