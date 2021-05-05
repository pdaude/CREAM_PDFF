function [A psif ampl freq] = fat_basis(te,Tesla,NDB,H2O,units)
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
if ~exist('NDB','var') || isempty(NDB)
    NDB = 2.5;
end
if ~exist('H2O','var') || isempty(H2O)
    H2O = 4.7;
end
if ~exist('units','var') || isempty(units)
    units = 'proton';
end

%% triglyderide properties

if NDB == -1
    % backward compatibility
    d    = [1.300 2.100 0.900 5.300 4.200 2.750];
    ampl = [0.700 0.120 0.088 0.047 0.039 0.006];
    awater = 1;
else
    % fat chemical shifts in ppm
    d = [5.29 5.19 4.20 2.75 2.24 2.02 1.60 1.30 0.90];
    
    % Mark's heuristic formulas
    CL = 16.8+0.25*NDB;
    NDDB = 0.093*NDB^2;
    
    % Gavin's formulas (no. protons per molecule)
    awater = 2;
    ampl(1) = NDB*2;
    ampl(2) = 1;
    ampl(3) = 4;
    ampl(4) = NDDB*2;
    ampl(5) = 6;
    ampl(6) = (NDB-NDDB)*4;
    ampl(7) = 6;
    ampl(8) = (CL-4)*6-NDB*8+NDDB*2;
    ampl(9) = 9;
    
    % the above counts no. molecules so that 1 unit of water = 2
    % protons and 1 unit of fat = sum(a) = (2+6*CL-2*NDB) protons.
    if isequal(units,'mass')
        % scale in terms of molecular mass so that 18 units of
        % water = 2 protons and (134+42*CL-2*NDB) units of fat =
        % (2+6*CL-2*NDB) protons. I.e. w and f are in mass units.
        awater = awater/18;
        ampl = ampl/(134+42*CL-2*NDB);
    else
        % scale by no. protons (important for PDFF) so 1 unit
        % of water = 1 proton and 1 unit of fat = 1 proton.
        awater = awater/2;
        ampl = ampl/sum(ampl);
    end
end

%% fat water matrix

% put variables in the right format
ne = numel(te);
te = reshape(te,ne,1);
te = gather(te);
Tesla = cast(Tesla,'like',te);

% time evolution matrix
fat = te*0;
water = te*0+awater;
larmor = 42.57747892*Tesla; % larmor freq (MHz)
for j = 1:numel(d)
    freq(j) = larmor*(d(j)-H2O); % relative to water
    fat = fat + ampl(j)*exp(2*pi*i*freq(j)*te);
end
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