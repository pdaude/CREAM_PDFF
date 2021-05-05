function h = ims(im,CLIM,TITLE)
%Show images: h = ims(im,CLIM,TITLE)
%
% CLIM:
%   [FRAC] = scale to fraction of the range
%   [CLOW CHIGH] = same as in imagesc
%   [0 1; NaN 100; ...] = individual ranges (NaN = no range)
%
% TITLE:
%   {'plot1', 'plot2', ...}
%
% h returns handles to all the subplots

% add flexibility to TITLE (e.g. can use 1:N as titles)
if exist('TITLE','var')
    if isnumeric(TITLE)
        TITLE = num2str(reshape(TITLE,[],1));
    else
        % ignore certain formatting characters
        TITLE = strrep(TITLE,'_','\_');
    end
    if iscell(TITLE)
        TITLE = reshape(TITLE,[],1);
    end
else
    TITLE = '';
end
if ~exist('CLIM','var') || isempty(CLIM)
    CLIM = 0.99; % default to 99% of range
end
if numel(CLIM)==1
    if CLIM<0 || CLIM>1 || ~isreal(CLIM)
        error('CLIM (scalar) must be between 0 and 1.');
    end
else
    if size(CLIM,2)~=2
        error('CLIM must have shape [low high; ... ].')
    end
end

% for complex produce 2 maps
if ~isreal(im)
    figure; ims(angle(im),[-pi pi],'phase'); % phase
    figure; im = abs(im); % magnitude
end

% make into 3D array
im = squeeze(im);
[x y n] = size(im);
im = reshape(im,x,y,n);

% catch likely errors
nmax = 8;
if n>nmax
    mid = round(n/2-nmax/2);
    warning('too many images (%i)... showing slices %i-%i only.',n,mid,mid+nmax)
    im = im(:,:,mid:mid+nmax);
    n = nmax;
end
if n==0
    error('Image array cannot be empty.')
end
if isa(im,'gpuArray')
    im = gather(im);
end

% try and figure out a nice tiling
rows = floor(sqrt(n));
cols = ceil(n/rows);
N = [rows cols];

% clear existing plots, unless a single image
if n>1
    clf reset
end

% plot here
for i = 1:N(1)
    for j = 1:N(2)
        k = (i-1)*N(2)+j;
        if k>n; break; end
        
        im_k = double(im(:,:,k)); % prevent any matlab weirdness
        lo = nanmin(nanmin(im_k));
        hi = nanmax(nanmax(im_k));
        
        if ~isempty(CLIM) && isfinite(lo) && isfinite(hi)
            if numel(CLIM)==1 % use fraction of the range
                edges = linspace(lo,hi,256);
                C = histc(reshape(im_k,[],1),edges);
                C = cumsum(C) / numel(im_k);
                [C IA] = unique(C);
                if C(1)==0
                    CP = C; IAP = IA;
                else
                    CP = [0;C]; IAP = [1;1+IA]; % pad a zero at the bottom edge
                end
                idx = interp1(CP,IAP,(1-CLIM)/2);
                lo = interp1(IA,edges(IA),idx);
                idx = interp1(CP,IAP,1-(1-CLIM)/2);
                hi = interp1(IA,edges(IA),idx);
            elseif size(CLIM,1)==1 % use same setting for all
                if ~isnan(CLIM(1,1));lo = CLIM(1,1);end
                if ~isnan(CLIM(1,2));hi = CLIM(1,2);end
            elseif k<=size(CLIM,1) % individual CLIMs
                if ~isnan(CLIM(k,1));lo = CLIM(k,1);end
                if ~isnan(CLIM(k,2));hi = CLIM(k,2);end
            end
        end
        % retile, unless a single image
        if n>1
            h(k) = subplot(N(1),N(2),k);
        end
        if lo<hi
            imagesc(im_k,[lo hi]);
        else
            imagesc(im_k);
        end
        set(gca,'XTickLabel','');
        set(gca,'YTickLabel','');
        if ~isempty(TITLE)
            if size(TITLE,1)==1 % use same setting for all
                title(TITLE(1,:),'FontSize',10)
            else
                if k<=size(TITLE,1) % individual titles
                    title(TITLE(k,:),'FontSize',10)
                end
            end
        else
            if exist('mid','var') % user sent too many images, label them
                title(mid+k-1,'FontSize',10)
            end
        end
    end
end

% if extra space, put focus on "next" plot
if k<=N(1)*N(2) && n>1
    h(k) = subplot(N(1),N(2),k);
    axis off
end

% prevent output to screen
if nargout==0; clear h; end
