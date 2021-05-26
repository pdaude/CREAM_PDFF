function [A psif ampl freq] = fat_basis_CREAM(te,Tesla,algoParams)
% [A psif ampl freq] = fat_basis(te,Tesla,NDB,H2O,units)
%
% Function that produces the fat-water matrix. Units can be adjusted
% to reflect the relative masses of water and fat. Default is proton.
% For backward compatibility, set NDB=-1 to use old Hamilton values.
%
% Inputs:
%  te (echo times in sec)
%  Tesla (field strength)
%  NDB (number of double bonds)
%  H2O (water freq in ppm)
%  units ('proton' or 'mass')
%
% Ouputs:
%  A = matrix [water fat] of basis vectors
%  psif = best fit of fat to exp(i*B0*te-R2*te)
%         where B0=real(psif) and R2=imag(psif)
%         note: B0 unit is rad/s (B0/2pi is Hz)
% ampl = vector of relative amplitudes of fat peaks
% freq = vector of frequencies of fat peaks (unit Hz)
%
% Ref: Bydder M, Girard O, Hamilton G. Magn Reson Imag. 2011;29:1041

%% argument checks

if max(te)<1e-3 || max(te)>1
    error('''te'' should be in seconds.');
end

%% fat water matrix
disp(Tesla)
% put variables in the right format
ne = numel(te);
te = reshape(te,ne,1);

larmor = algoParams.gyro*Tesla; % larmor freq (MHz)

H2O=algoParams.species(1).frequency; % 4.7
water=algoParams.species(1).relAmps*ones(ne,1);
ampl=algoParams.species(2).relAmps;

freq=larmor*(algoParams.species(2).frequency-algoParams.species(1).frequency);
FatA = ampl.*exp(2*pi*i*te*freq);
fat=sum(FatA,2);

A = [water fat];

%% nonlinear fit of fat to complex exp (gauss newton)

if nargout>1
    psif = [2*pi*larmor*(1.3-H2O) 50]; % initial estimates (rad/s)
    psif = double(gather(psif)); % keep MATLAB happy
    psif = fminsearch(@(psif)myfun(psif,double(te),double(fat)),psif);
    psif = cast(complex(psif(1),psif(2)),'like',te);
    %te2=linspace(0,max(te),100*numel(te)); cplot(1000*te,A(:,2),'o');
    %hold on; cplot(1000*te2,exp(i*psif*te2)); title(num2str(psif)); hold off; keyboard
end

% exponential fitting function
function normr = myfun(psif,te,data)
psif = complex(psif(1),psif(2));
f = exp(i*psif*te); % function
v = (f'*data)/(f'*f); % varpro
r = v*f-data; % residual
normr = double(gather(norm(r)));