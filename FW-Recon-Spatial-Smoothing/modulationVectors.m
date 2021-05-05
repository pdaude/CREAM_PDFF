function [B, Bh] = modulationVectors(nB0, N)
%MODULATIONVECTORS Returns a matrix used for modelling the signal
%   Returns a matrix used for modelling the signal at different
%   off-resonance frequencies

B = zeros(nB0,N,N);
Bh = zeros(nB0,N,N);
for b = 0:nB0-1
    omega = 2*pi*b/nB0;
    B(b+1,:,:) = eye(N);
    for n = 0:N-1
        B(b+1,n+1,n+1) = exp(1i*n*omega);
    end
    Bh(b+1,:,:) = conj(B(b+1,:,:));
end

end