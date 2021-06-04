function [sse,YM]=calculate_residual(params,algoParams,imDataParams)
%params.B0 % B0 (Hz)
%params.R2 % R2* (1/s)
%params.FF % FF (%)
%params.PH % initial phase (rad)
%params.F  Fat image;
%params.W Water image;
R2=reshape(params.R2,[1,numel(params.R2)]);% [1;nbrVoxels]
B0=reshape(params.B0,[1,numel(params.B0)]);% [1;nbrVoxels]
F=reshape(params.F,[1,numel(params.F)]);% [1;nbrVoxels]
W=reshape(params.W,[1,numel(params.W)]);% [1;nbrVoxels]
TE=imDataParams.TE; % second
PHI=zeros(1,numel(params.B0));
if isfield(params,'PH')
   PHI= reshape(params.PH,[1,numel(params.PH)]);
end
[nx,ny,nz,ncoils,ne]=size(imDataParams.images);
if ncoils >1 
        error('''varargin'' must be option/value pairs.');
end
image = reshape(imDataParams.images,nx,ny,nz,ne);
% echos in 1st dim
image = permute(image,[4 1 2 3]);
image = reshape(image,ne,nx*ny*nz);

%% fat water amplitudes matrix 
larmor = algoParams.gyro*imDataParams.FieldStrength; % larmor freq (MHz)
FreqA=(algoParams.species(2).frequency-algoParams.species(1).frequency);
FatA = algoParams.species(2).relAmps.*exp(2*pi*i*larmor*TE'*FreqA);

A = [ones(numel(TE),1), sum(FatA,2)];



%YM=B*A*X
% B= exp(i(B0te+PHI)+R2te)
B = exp(i*(TE'*complex(2*pi*B0,R2)+repmat(PHI,numel(TE),1)));
X=[W;F];

YM=B.*(A*X);
r = reshape(image-YM,ne,nx,ny,nz);
sse = squeeze(real(dot(r,r))); % sum of squares error
YM=reshape(YM,ne,nx,ny,nz);
YM = permute(YM,[2 3 4 1]);
end