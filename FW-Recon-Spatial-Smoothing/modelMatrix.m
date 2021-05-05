function RA = modelMatrix(N, t1, dt, M, R2, gyro, B0, alpha, CS, P)
%MODELMATRIX Returns a matrix used for modelling the signal
%   Returns a matrix used for modelling the signal, given specified
%   echo times, R2*, and field strength
    RA = zeros(N, M);
    for n = 0:N-1
        t = t1+n*dt;
        RA(n+1, 1) = exp(-(t-t1)*R2);  % Water resonance
        for p = 1:P-1  % Loop over fat resonances
            % Chemical shift between water and the different fat peaks (in ppm)
            omega = 2*pi*gyro*B0*(CS(p+1)-CS(1));
            RA(n+1, 2) = RA(n+1, 2) + alpha(2,p+1)*exp(-(t-t1)*R2 + 1i*t*omega);
        end
    end
end