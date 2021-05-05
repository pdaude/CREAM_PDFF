function [params sse] = ndbest(te,data,Tesla,H2O)
%[params sse] = ndbest(te,data,Tesla,H2O)
%
% Fat water separation with double bond (NDB) estimation.
% WARNING - VERY SLOW!! Only does single pixel fitting.
%
% Inputs:
%  te is a vector in seconds (ne)
%  data is a complex vector (ne)
%  Tesla is field strength (scalar)
%  H2O is the water ppm (optional)
%
% Outputs:
%  params.B0 is B0 (Hz)
%  params.R2 is R2* (1/s)
%  params.FF is FF (%)
%  params.PH is PH (rad)
%  params.NDB is NDB

display = 1; % slow but useful

%% sizes

ne = numel(te);
te = reshape(te,ne,1);
data = reshape(data,ne,1);

%% initialize parameters

% initial estimates
tmp = angle(dot(data(1:ne-1),data(2:ne)));
B0 = tmp/mean(diff(te)); % B0 in rad/s
R2 = 50; % R2* in 1/s
NDB = 3; % no. double bonds

% bounds (B0 R2 NDB)
LB = [-Inf  0  1];
UB = [+Inf 200 6];

% water freq in ppm
if nargin<4 || isempty(H2O)
    H2O = 4.7; 
end

% make lsqnonlin happy
te = double(gather(te));
data = double(gather(data));
Tesla = double(gather(Tesla));
B0 = double(gather(B0));
R2 = double(gather(R2));
NDB = double(gather(NDB));
H2O = double(gather(H2O));

%% nonlinear fitting

opts = optimset('display','off');

% first estimate (assume B0 is on water peak)
x1 = [B0 R2 NDB];
[x1,sse1,~,~,~,~,J1] = lsqnonlin(@(x)myfun(x,te,data,Tesla,H2O),x1,LB,UB,opts);

% second estimate (assume B0 is on fat peak)
[~,psif] = fat_basis(te,Tesla,NDB,H2O);
x2(1) = x1(1)-real(psif);
x2(2) = max(x1(2)-imag(psif),0);
x2(3) = NDB;
[x2,sse2,~,~,~,~,J2] = lsqnonlin(@(x)myfun(x,te,data,Tesla,H2O),x2,LB,UB,opts);

% choose best solution
if sse1 < sse2
    x = x1;
    J = J1;
    sse = sse1;
else
    x = x2;
    J = J2;
    sse = sse2;
end

% get the linear(ish) terms (wf p)
[r wf] = myfun(x,te,data,Tesla,H2O);

% display
if display
    
    cplot(1000*te,data,'o'); xlabel('te (ms)');
    te2 = linspace(0,min(te)+max(te),10*ne);
    data2 = myfun([x(:);wf(:)],te2,data,Tesla,H2O);
    hold on; cplot(1000*te2,data2); hold off
    title(['FF=' num2str(100*real(wf(2)/sum(wf)),'%.1f') ' NDB=' num2str(x(3),'%.2f')]);
    axis tight; drawnow

end

% return variable
params.B0 = x(1)/2/pi; % units Hz
params.R2 = x(2); % units 1/s
params.FF = min(max(100*real(wf(2)/sum(wf)),-10),110); 
params.PH = angle(sum(wf));
params.NDB = x(3);

% display
if display
    
    v = 2*ne-6; % no. degrees of freedom
    cov = pinv(full(J'*J))*sse/v; % covariance matrix
    ci95 = sqrt(max(diag(cov),0))*1.96; % confidence intervals
    
    disp([' initial B0 ' num2str(B0/2/pi) ' R2* ' num2str(R2)])
    disp([' B0    ' num2str(x(1)/2/pi) ' ± ' num2str(ci95(1)/2/pi)])
    disp([' R2*   ' num2str(x(2)) ' ± ' num2str(ci95(2))])
    disp([' FF    ' num2str(params.FF)])
    disp([' PH    ' num2str(params.PH)])
    disp([' NDB   ' num2str(x(3)) ' ± ' num2str(ci95(3))])
    disp([' H2O   ' num2str(H2O) ' (fixed)'])
    disp([' sse   ' num2str(sse,10)])
    disp(' ')
    
end

%% function to calculate residual
function [r wf] = myfun(x,te,data,Tesla,H2O)

% unknowns
B0 = x(1);  % field map (rad/s)
R2 = x(2);  % exp decay (1/s)
NDB = x(3); % no. double bonds

% water/fat time evolution matrix
A = fat_basis(te,Tesla,NDB,H2O);

% fieldmap and R2*
W = diag(exp(i*complex(B0,R2)*te));

% two paths: one for residual, one for display
if numel(x)==3

    % calculate water, fat and initial phase
    Mh = pinv(real(A'*W'*W*A));
    Ab = A'*W'*data;
    p = angle(sum(Ab.*(Mh*Ab)))/2;
    wf = Mh*real(Ab*exp(-i*p))*exp(i*p);

    % recalculate p to absorb sign of wf
    if nargout>1
        p = angle(sum(wf));
        wf = real(wf*exp(-i*p))*exp(i*p);
    end
    
    % residual
    r = reshape(W*A*wf-data,size(te));
    r = [real(r);imag(r)];

else

    % display: function values from provided x
    wf = x(4:5);
    r = W*A*wf;

end
