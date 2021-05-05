function [water, fat, FF, R2map, B0map ,label, lambda_map] = doReconGauss( load_file, voxelsize, param )
%DO_RECON_GAUSS Main file for reconstruction
%   Returns water image
%   Returns fat image
%   Returns FF image, Unit: per mille
%   Returns R2map, Unit: s^-1
%   Returns B0map, Unit: ppm
%   Returns label image, tells which label was assigned by QPBO, for debugging
%   Returns lambda_map, tells which lambda value was used in the final iteration, for debugging

nR2 = param.nR2;
R2max = param.R2max;
R2cand = param.R2cand;
lambda = param.lambda;
nB0 = param.nB0;
nICMiter = param.nICMiter;
use3D = param.use3D;
maxICMupdate = param.maxICMupdate;
sigma = param.sigma;

fatCS = param.fatCS;
relAmps = param.relAmps;
watCS = param.watCS;

gyro = param.gyro;

addpath(['.' filesep 'qpboMex-master'])

R2step = R2max/(nR2-1);

M = 2;%number of species (water and fat)
CS = [watCS fatCS];
P = length(CS);

relAmps = relAmps/sum(relAmps);%normalization

alpha = zeros(M,P);
alpha(1,1) = 1;
alpha(2,2:end) = relAmps;

dx = voxelsize(1);
dy = voxelsize(2);
dz = voxelsize(3);

load(load_file)

if imDataParams.PrecessionIsClockwise ~= 1
    imDataParams.images = real(imDataParams.images) - 1i*imag(imDataParams.images);
end

TE = imDataParams.TE;
B0 = imDataParams.FieldStrength;

t1 = TE(1);
dt = mean(diff(TE));

Y = imDataParams.images;
sz = size(Y);
nx = sz(2);
ny = sz(1);
nz = sz(3);
Y = squeeze(Y);
Y = permute(Y,[4 2 1 3]);

Y = Y(:,:);
sz = size(Y);

nVxl = sz(2);
N = sz(1);
[B, Bh] = modulationVectors(nB0, N);

RA = zeros(nR2,N,M);
RAp = zeros(nR2,M,N);
for r = 0:nR2-1
    R2 = r*R2step;
    RA(r+1,:,:) = modelMatrix(N, t1, dt, M, R2, gyro, B0, alpha, CS, P);
    RAp(r+1,:,:) = pinv(squeeze(RA(r+1,:,:)));
end

C = zeros(nR2,nB0,N,N);
Qp = zeros(nR2,nB0,M,N);
for r = 0:nR2-1
    % Null space projection matrix
    proj = eye(N) - squeeze(RA(r+1,:,:)) * squeeze(RAp(r+1,:,:));
    for b = 0:nB0-1
        C(r+1,b+1,:,:) = squeeze(B(b+1,:,:)) * proj * squeeze(Bh(b+1,:,:));
        Qp(r+1,b+1,:,:) = squeeze(RAp(r+1,:,:)) * squeeze(Bh(b+1,:,:));
    end
end

% For B0 index -> off-resonance in ppm
B0step = 1.0/nB0/dt/gyro/B0;
V = zeros(nB0,1);  % Precalculate discontinuity costs
for b = 0:nB0-1
    V(b+1) = min(b^2, (b-nB0)^2);
end

J = getB0Residuals(Y, B, Bh, nB0, nVxl, R2cand, N, M, t1, dt, gyro, B0, alpha, CS, P);
J_orig = J;

if sigma
    H = getKernel(dx, dy, dz, sigma, use3D);
    grej = ones(nx,ny,nz);
    grej = imfilter(grej,H);
    for i = 1:nB0
        qq = reshape(J(i,:),[nx ny nz]);
        qq = imfilter(qq,H)./grej;
        J(i,:) = qq(:);
    end
end

[dB0, label, lambda_map] = calculateFieldMap(nB0, J, V, lambda, dx, dy, dz, nx, ny, nz, use3D, maxICMupdate, nICMiter, J_orig);

J = getR2Residuals(Y, dB0, C, nB0, nR2, nVxl);
R2 = greedyR2(J, nVxl);

% Find least squares solution given dB0 and R2
rho = zeros(M, nVxl);
for r = 0:nR2-1
    for b = 1:nB0
        vxls = (dB0 == b) & (R2 == r);
        y = Y(:, vxls);
        rho(:, vxls) = squeeze(Qp(r+1,b,:,:)) * y;
    end
end

R2map = R2*R2step;
B0map = (dB0-1)*B0step;

R2map = permute(reshape(R2map,[nx ny nz]),[2 1 3]);
B0map = permute(reshape(B0map,[nx ny nz]),[2 1 3]);
water = permute(reshape(rho(1,:),[nx ny nz]),[2 1 3]);
fat = permute(reshape(rho(2,:),[nx ny nz]),[2 1 3]);

water = abs(water);
fat = abs(fat);

FF = 1000*(fat)./(fat + water);
FF(isnan(FF)) = 0;

label = permute(reshape(label,[nx ny nz]),[2 1 3]);
lambda_map = permute(reshape(lambda_map,[nx ny nz]),[2 1 3]);

end