function B = medfiltn(A,P,mask)
%
% Median filter of N-D array using kernel P.
% Performs circular wrapping at boundaries.
%
% A is an N-D array
% P is a kernel (default [3 3])
% mask (binary) excludes points

% argument checks
if nargin<2 || isempty(P)
    P = [3 3];
elseif ~isvector(P) || any(P<1) || nnz(mod(P,1))
    error('something wrong with P');
end
if ~isreal(A)
    error('median not well defined for complex.');
end
if numel(P)>ndims(A)
    error('P has more dims than A');
else
    P = reshape(P,numel(P),1);
end
if nargin>2
    if ~isequal(size(A),size(mask)) || nnz(mod(mask,1))
        error('something wrong with mask');
    end
    origA = A;
    A(~mask) = NaN; % let's use use omitnan flag
end

% generate shift indicies
S = zeros(numel(P),prod(P));

for k = 1:prod(P)
    switch numel(P)
        case 1; [S(1,k)] = ind2sub(P,k);
        case 2; [S(1,k) S(2,k)] = ind2sub(P,k);
        case 3; [S(1,k) S(2,k) S(3,k)] = ind2sub(P,k);
        case 4; [S(1,k) S(2,k) S(3,k) S(4,k)] = ind2sub(P,k);
        otherwise; error('high dim not implemented - fix me');
    end
end

% make data matrix
B = zeros(numel(A),prod(P),'like',A);

for k = 1:prod(P)
    tmp = circshift(A,S(:,k)-fix(P/2)-1);
    B(:,k) = reshape(tmp,numel(A),1);
end

% median along shift dim
B = median(B,2,'omitnan');

% use orig values when all masked
if nargin>2
    k = isnan(B);
    B(k) = origA(k);
end

% return it as we found it
B = reshape(B,size(A));
