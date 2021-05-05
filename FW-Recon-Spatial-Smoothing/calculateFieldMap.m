function [dB0, label, lambda_map] = calculateFieldMap(nB0, J, V, lambda, dx, dy, dz, nx, ny, nz, use3D, maxICMupdate, nICMiter, J_orig)
%CALCULATEFIELDMAP Calculates field map
%   Calculates field map

[A, B] = findTwoSmallestMinima(J);
dB0 = A;

% Prepare MRF
% Prepare discontinuity costs
vxls = size(J,2);
% 2nd derivative of residual function
% NOTE: No division by square(steplength) since
% square(steplength) not included in V
ddJ = J(sub2ind(size(J),mod(A, nB0)+1,(1:vxls)')) + J(sub2ind(size(J),mod(A-2, nB0)+1,(1:vxls)')) - 2*J(sub2ind(size(J),A,(1:vxls)'));

[left, right, up, down, above, below] = getIndexImages(nx, ny, nz);

sz1 = (nx-1) * ny * nz;
sz2 = nx * (ny-1) * nz;
sz3 = nx * ny * (nz-1);
edgeWeights = zeros(sz1 + sz2 + sz3,6);

ind = reshape(1:vxls,[nx ny nz]);

% left-right
ind1 = ind(left);
ind2 = ind(right);
wx = min(ddJ(left),ddJ(right))/dx;
w1 = V(abs(A(left)-A(right))+1);
w2 = V(abs(A(left)-B(right))+1);
w3 = V(abs(B(left)-A(right))+1);
w4 = V(abs(B(left)-B(right))+1);
edgeWeights(1:sz1,:) = [ind1(:) ind2(:) repmat(wx,[1 4]).*[w1(:) w2(:) w3(:) w4(:)]];

% up-down
ind1 = ind(up);
ind2 = ind(down);
wy = min(ddJ(up),ddJ(down))/dy;
w1 = V(abs(A(up)-A(down))+1);
w2 = V(abs(A(up)-B(down))+1);
w3 = V(abs(B(up)-A(down))+1);
w4 = V(abs(B(up)-B(down))+1);
edgeWeights(sz1+1:sz1+sz2,:) = [ind1(:) ind2(:) repmat(wy,[1 4]).*[w1(:) w2(:) w3(:) w4(:)]];

% above-below
if use3D
    ind1 = ind(above);
    ind2 = ind(below);
    wz = min(ddJ(above),ddJ(below))/dz;
    w1 = V(abs(A(above)-A(below))+1);
    w2 = V(abs(A(above)-B(below))+1);
    w3 = V(abs(B(above)-A(below))+1);
    w4 = V(abs(B(above)-B(below))+1);
    edgeWeights(sz1+sz2+1:sz1+sz2+sz3,:) = [ind1(:) ind2(:) repmat(wz,[1 4]).*[w1(:) w2(:) w3(:) w4(:)]];
else
    edgeWeights = edgeWeights(1:sz1+sz2,:);
end

% Prepare data fidelity costs
terminalWeights = lambda * [J(sub2ind(size(J),A,(1:vxls)')),J(sub2ind(size(J),B,(1:vxls)'))];

% QPBO
try
    [~, label] = qpboMex(terminalWeights, edgeWeights);
catch ME
    if (strcmp(ME.identifier,'MATLAB:scriptNotAFunction'))
        msg = 'Did not find compiled version of QPBO for your OS.';
        causeException = MException('QPBO:notCompiled',msg);
        ME = addCause(ME,causeException);
    end
    rethrow(ME)
end

% Method to fix unlabeled voxels
lambda_map = zeros(nx, ny, nz) + lambda;
% Increasing lambda
while find(label(:) == -1)
    terminalWeights(label==-1,:) = 2*terminalWeights(label==-1,:);
    lambda_map(label==-1) = 2*lambda_map(label==-1);
    [~, label] = qpboMex(terminalWeights, edgeWeights);
end

dB0(label == 0) = A(label == 0);
dB0(label == 1) = B(label == 1);

if nICMiter
    if ~use3D
        wz = [];
    end
    J = J_orig;
    [A, ~] = findTwoSmallestMinima(J);
    ddJ = J(sub2ind(size(J),mod(A, nB0)+1,(1:vxls)')) + J(sub2ind(size(J),mod(A-2, nB0)+1,(1:vxls)')) - 2*J(sub2ind(size(J),A,(1:vxls)'));
    wx = min(ddJ(left),ddJ(right))/dx;
    wy = min(ddJ(up),ddJ(down))/dy;
    if use3D
        wz = min(ddJ(above),ddJ(below))/dz;
    end
    dB0 = ICM(dB0, nB0, maxICMupdate, nICMiter, J, V, wx, wy, wz, left, right, up, down, above, below, lambda);
end

end