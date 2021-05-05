function [ fieldmapFine ] = findRangeMin( imDataParams, algoParams,fieldmapRough, r2starmap )
images = imDataParams.images;

try
  precessionIsClockwise = imDataParams.PrecessionIsClockwise;
catch
  precessionIsClockwise = 1;
end

  
% If precession is clockwise (positive fat frequency) simply conjugate data
if precessionIsClockwise <= 0 
  imDataParams.images = conj(imDataParams.images);
  imDataParams.PrecessionIsClockwise = 1;
end

gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:))*(imDataParams.FieldStrength)];
relAmps = algoParams.species(2).relAmps;
%range_fm = algoParams.range_fm;
t = imDataParams.TE;
%NUM_FMS = algoParams.NUM_FMS; 
Gridsize = algoParams.gridsize;
range_r2star = algoParams.range_r2star;
NUM_R2STARS = algoParams.NUM_R2STARS;


% Images are size sx X sy, N echoes, C coils
sx = size(images,1);
sy = size(images,2);
N = size(images,5);
C = size(images,4);

% Get VARPRO-formulation matrices for given echo times and chemical shifts 
Phi = getPhiMatrixMultipeak(deltaF,relAmps,t);
%Phi(:,1) = 0;
iPhi = pinv(Phi'*Phi);
A = Phi*iPhi*Phi';

fieldmapRough = round(fieldmapRough);
tempfield = fieldmapRough(:);
%maxfield = max(tempfield);
%minfield = min(tempfield);
temp = reshape(squeeze(images), [sx*sy, N]);
temp = permute(temp, [2,1]); % 6by192*192

for k=1:length(tempfield)
    P1 = [];
    searchrange = tempfield(k)-round(Gridsize/2) : tempfield(k)+round(Gridsize/2);
    for s = 1:length(searchrange)
        Psi = diag(exp(1i*2*pi*searchrange(s)*t - abs(t)*r2starmap(k)));
        mm = (eye(N)-Psi*Phi*pinv(Psi*Phi));
        P1 = [P1; mm];
    end
    res = sum(abs(reshape(P1*temp(:,k), [N, length(searchrange)])).^2,1);
    [~,ind] = min(res,[], 2);
    fieldmapFine(k) = searchrange(ind);
end
fieldmapFine = reshape(fieldmapFine,[sx, sy]);

% Compute residual
% r2s =linspace(range_r2star(1),range_r2star(2),NUM_R2STARS);
% 
% % Precompute all projector matrices (one per field value) for VARPRO
% fiedmapRange = minfield-Gridsize:maxfield+Gridsize;
% P = [];
% for kr=1:NUM_R2STARS
%     P1 = [];
%   for k=1:length(fiedmapRange)
%     Psi = diag(exp(1i*2*pi*fiedmapRange(k)*t - abs(t)*r2s(kr)));
%     mm = (eye(N)-Psi*Phi*pinv(Psi*Phi));
%     P1 = [P1; mm];
% 
%   end
%     P(:,:,kr) = P1;  %  4260by6by11
% end
% P = permute(P,[3,1,2]); %11by 4260by6
% 
% 
% temp = reshape(squeeze(images), [sx*sy, N]);
% temp = permute(temp, [2,1]); % 6by192*192
% temp2 = zeros(1,length(tempfield));
% 
% for  kr=1:NUM_R2STARS
%     residue_temp1(kr,:,:) = squeeze(P(kr,:,:))*temp;
% end
% % residue 11by 4260by 192*192
% for t = 1:length(tempfield)
%     ind_r2starmap = find(r2starmap(t)==r2s);
%     residue_temp2(:,t) = squeeze(residue_temp1(ind_r2starmap,:,t)); %4260by 192*192
% %     residue = squeeze(P(ind_r2starmap,:,:))*temp(:,t); %4260by1
% %     residue = reshape(squeeze(residue),[N,length(fiedmapRange)]);% 6byRange
% %     residue = sum(abs(residue).^2,1);
% %     [nnn,ind] = min(residue( tempfield(t)-Gridsize+abs(minfield-Gridsize)+1:tempfield(t)+Gridsize+abs(minfield-Gridsize)+1  ));
% %     %residue = temp(t,:)*P(:,(tempfield(t)-Gridsize):(tempfield(t)+Gridsize),r2starmap(t)); %1byN times N by 2*Gridsize+1; N is the echoes of time;
% %     %[~,ind] = min(residue);
% %     temp2(t) = tempfield(t)+ind - (Gridsize+1);
% end
% residue = reshape(residue_temp2,[N,length(fiedmapRange), sx*sy]);% 6byRangeby sx*sy
% Res = squeeze(sum(residue, 1).^2);
% % search range for each pixel in the fieldmap
% searchrange = round(Gridsize/2);
% lowbound = -searchrange+abs(minfield-Gridsize)+1;
% upbound = searchrange+abs(minfield-Gridsize)+1;
% for t = 1:length(tempfield)
%     [~,ind] = min(Res(tempfield(t)+lowbound:tempfield(t)+upbound,t),1);
%     temp2(t) = tempfield(t)+ind - (searchrange+1);
% end
% 
% fieldmapFine = reshape(temp2,[sx, sy]);

end

