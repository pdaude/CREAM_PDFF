function [params sse] = lipoquant(te,data,Tesla,varargin)
%[params sse] = lipoquant(te,data,Tesla,varargin)
%
% Proton density fat fraction estimation using magnitude.
% Returns struct of FF (%) and R2* (1/s).
%
% Inputs:
%  te is echo time in seconds (vector)
%  data (1D, 2D or 3D with echos in last dimension)
%  Tesla field strength (scalar)
%  varargin in option/value pairs (e.g. 'init',0.1)
%
% Alternate input:
%  imDataParams is a struct from fat-water toolbox
%
% Outputs:
%  params.R2 is R2* (1/s)
%  params.FF is FF (%)
%  sse is sum of squares error

% demo code
if nargin==0
    load phantom_16echo_bipolar.mat
    %load PHANTOM_NDB_PAPER.mat
    %load liver_gre_3d_2x3.mat; data = data(:,:,31,:);
end

%% default options

opts.ndb = 2.5; % no. double bonds
opts.h2o = 4.7; % water ppm
opts.maxit = 20; % max. no. iterations
opts.nonnegR2 = 1; % nonnegative R2*
opts.noise = 0; % noise std (if available)
opts.init = 0.1; % initial FF estimate (0-1)

% varargin handling (must be option/value pairs)
for k = 1:2:numel(varargin)
    if k==numel(varargin) || ~ischar(varargin{k})
        error('''varargin'' must be option/value pairs.');
    end
    if ~isfield(opts,varargin{k})
        error('''%s'' is not a valid option.',varargin{k});
    end
    opts.(varargin{k}) = varargin{k+1};
end

%% special handling to accept ISMRM F/W toolbox structure

if isa(te,'struct')
    Tesla = te.FieldStrength;
    [nx ny nz nc ne] = size(te.images);
    if nc>1; error('multiple coils not supported'); end
    data = reshape(te.images,nx,ny,nz,ne);
    te = te.TE; % seconds
    data = abs(data);
end

%% argument checks - be flexible: 1D, 2D or 3D data

if numel(data)==numel(te)
    data = reshape(data,1,1,1,[]);
elseif ndims(data)==2 && size(data,2)==numel(te)
    data = permute(data,[1 3 4 2]);
elseif ndims(data)==3 && size(data,3)==numel(te)
    data = permute(data,[1 2 4 3]);
elseif ndims(data)==4 && size(data,4)==numel(te)
else
    error('data [%s\b] not compatible with te [%s\b].',...
    sprintf('%i ',size(data)),sprintf('%i ',size(te)));
end
if numel(te)<3
    warning('too few echos (%i) for FF, R2* estimation.',numel(te));
end
if max(te)>1
    error('''te'' should be in seconds.');
end
if ~isfloat(data)
    warning('data should be floating point. Converting to double.');
    data = double(data);
end
if ~isreal(data) || nnz(data<0)
    warning('data should be magnitude. Converting to magnitude.');
    data = abs(data);
end
[nx ny nz ne] = size(data);

%% see if gpu is possible

try
    gpu = gpuDevice;
    if verLessThan('matlab','8.4'); error('GPU needs MATLAB R2014b'); end
    data = gpuArray(data);
    fprintf(' GPU found = %s (%.1f Gb)\n',gpu.Name,gpu.AvailableMemory/1e9);
catch ME
    data = gather(data);
    warning('%s. Using CPU.',ME.message);
end

% display
disp([' Data size: [' sprintf('%i ',size(data)) sprintf('\b]')])
disp([' TE (ms): ' sprintf('%.2f ',1000*te(:))])
disp([' Tesla (T): ' sprintf('%.2f',Tesla)])
disp(opts);

%% faster calculations

% echos in 1st dimension
data = permute(data,[4 1 2 3]);

% consistent data types
te = reshape(real(cast(te,'like',data)),ne,1);

%% nonlinear estimation (local optmization)

% initial estimates
if ~isreal(opts.init) || nnz(opts.init<0 | opts.init>1)
    warning('opts.init values should be [0 1].');
end
opts.init = min(max(real(opts.init),0),1);

if isscalar(opts.init)
    F = repmat(cast(opts.init,'like',data),1,nx*ny*nz);
else
    F = reshape(cast(opts.init,'like',data),1,nx*ny*nz);
end
M = reshape(max(data),1,nx*ny*nz);
R = repmat(cast(0,'like',data),1,nx*ny*nz);

% find the local minimum
[sse M F R] = nlsfit(M,F,R,te,data,Tesla,opts);

%% return parameters in correct format

% nonegative M and R2*
M = abs(M);
if opts.nonnegR2; R = abs(R); end

M = reshape(M,nx,ny,nz);
F = reshape(F,nx,ny,nz);
R = reshape(R,nx,ny,nz);
sse = reshape(sse,nx,ny,nz);

params.FF =  gather(100*F); % FF (%)
params.R2 = gather(R); % R2* (1/s)

sse = gather(sse);
fprintf(' sse = %.6e\n',sqrt(sum(sse(:))));

%% nonlinear least squares fitting
function [sse M F R] = nlsfit(M,F,R,te,data,Tesla,opts)

[ne nx ny nz] = size(data);

% water-fat matrix
opts.A = fat_basis(te,Tesla,opts.ndb,opts.h2o);

% residual and Jacobian
[r JM JF JR] = myfun(M,F,R,te,data,opts);

for iter = 1:opts.maxit

    % Normal matrix for the Jacobian is 3x3, which has an exact inverse:
    %
    % | a11 a12 a13 |-1           |   a33.a22-a32.a23 -(a33.a12-a32.a13)  a23.a12-a22.a13  |
    % | a21 a22 a23 |   = 1/dtm * | -(a33.a21-a31.a23)  a33.a11-a31.a13 -(a23.a11-a21.a13) |
    % | a31 a32 a33 |             |   a32.a21-a31.a22 -(a32.a11-a31.a12)  a22.a11-a21.a12  |
    %
    % with dtm = a11.(a33.a22-a32.a23)-a21.(a33.a12-a32.a13)+a31.(a23.a12-a22.a13)

    a11 = dot(JM,JM);
    a12 = dot(JM,JF);
    a13 = dot(JM,JR);
    a21 = a12;
    a22 = dot(JF,JF);
    a23 = dot(JF,JR);
    a31 = a13;
    a32 = a23;
    a33 = dot(JR,JR);

    % Levenberg Marquardt term
    if iter==1; lambda = (a11+a22+a33+mean(data(:)))/1000; end
 
    a11 = a11 + lambda;
    a22 = a22 + lambda;
    a33 = a33 + lambda;

    dtm = a11.*(a33.*a22-a32.*a23)-a21.*(a33.*a12-a32.*a13)+a31.*(a23.*a12-a22.*a13);
    
    % gradient (J'*r)
    gM = dot(JM,r);
    gF = dot(JF,r);
    gR = dot(JR,r);

    % updates = -inv(J'J+lambda.I)*J'r
    dM =(-(a33.*a22-a32.*a23).*gM+(a33.*a12-a32.*a13).*gF-(a23.*a12-a22.*a13).*gR)./dtm; 
    dF =(+(a33.*a21-a31.*a23).*gM-(a33.*a11-a31.*a13).*gF+(a23.*a11-a21.*a13).*gR)./dtm;
    dR =(-(a32.*a21-a31.*a22).*gM+(a32.*a11-a31.*a12).*gF-(a22.*a11-a21.*a12).*gR)./dtm;
    
    % evaluate at new estimate
    [r_new JM_new JF_new JR_new] = myfun(M+dM,F+dF,R+dR,te,data,opts);

    % accept changes that gave an improvement
    ok = dot(r_new,r_new) < dot(r,r);
    r(:,ok) = r_new(:,ok);
    M(ok) = M(ok)+dM(ok);
    F(ok) = F(ok)+dF(ok);
    R(ok) = R(ok)+dR(ok);
    JM(:,ok) = JM_new(:,ok);
    JF(:,ok) = JF_new(:,ok);
    JR(:,ok) = JR_new(:,ok);
   
    % update damping term
    lambda(ok) = lambda(ok) / 10;
    lambda(~ok) = lambda(~ok) * 9;

end

% error norm
sse = sum(r.^2);

% display maps
if numel(data)==ne
    plot(1e3*te,data,'o'); temp = linspace(min(te),max(te),100*numel(te)).';
    opts.A = fat_basis(temp,Tesla,opts.ndb,opts.h2o); xlabel('te (ms)');
    hold on; plot(1e3*temp,myfun(M,F,R,temp,data,opts)); hold off
    if opts.nonnegR2; tmp = abs(R); else tmp = R; end
    title(['FF=' num2str(100*F,'%.1f%%') '  R2=' num2str(tmp,'%.1f') '  ||r||=' num2str(sqrt(sse))]);
    fprintf('R2=%.3f\tFF=%.3f\tsse=%.5e\n',tmp,F,sse);pause(0.1);
else
    mid = ceil(nz/2);
    subplot(1,2,1); tmp = reshape(F,nx,ny,nz);
    imagesc(100*tmp(:,:,mid),[-5 55]); colorbar; title('PDFF (%)');
    subplot(1,2,2); tmp = reshape(R,nx,ny,nz);
    if opts.nonnegR2; tmp = abs(tmp); end
    imagesc(tmp(:,:,mid),[0 155]); colorbar; title('R2* (s^{-1})'); drawnow
end

%% model function
function [f JM JF JR] = myfun(M,F,R,te,data,opts)
  
% no. echos
ne = numel(te);

% nonnegative M and R2*
M = abs(M);
if opts.nonnegR2; R = abs(R); end

% components of the signal
fat = opts.A(:,2)*F;
water = opts.A(:,1)*[1-F];

s1 = repmat(M,ne,1);
s2 = abs(water+fat);
s3 = exp(-te*R);

% signal
f = hypot(s1.*s2.*s3,opts.noise);

% derivatives
if nargout>1
   
    JM = s1.*s2.^2.*s3.^2./f;
    JR = bsxfun(@times,-te,JM.*s1);
    
    if opts.nonnegR2
        signR = sign(R)+(R==0);
        JR = JR.*signR;
    end

    % JF (analytical)    
    tmp = abs(opts.A(:,2)).^2-real(opts.A(:,1)).^2-2*real(prod(opts.A,2));
    tmp = bsxfun(@plus,diff(real(opts.A),1,2),tmp*F);
    JF = s1.^2.*s3.^2.*tmp./f;  

    % JF (finite difference)
    %h = 1e-4;
    %s4 = bsxfun(@times,diff(opts.A,1,2),h);
    %s2 = abs(water+fat+s4);
    %JF = (hypot(s1.*s2.*s3,opts.noise)-f)./h;

end

% residual (or display)
if numel(f)==numel(data)
    f = f-reshape(data,size(f));
end
