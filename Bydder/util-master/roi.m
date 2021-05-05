function ind = roi(im,N)
%function roi()
% Interatively draw ROI on current figure (or on supplied image im).
% Displays number, mean and 95% confidence interval.
% Returns a cell array of vectors of the indices for each roi.
% If present, N is the number of ROIs to draw.

% user sets current figure
if nargin<1
    waitforbuttonpress;
else
    imagesc(im);
    title('Draw ROI')
end
if nargin<2
    N = Inf;
end

% make it big
h = gcf;
siz = get(h,'Position');
set(h,'Position', [100, 300, 900, 700]);

% loop until user selects an empty ROI (or N is reached)
ind = {};
while numel(ind)<N

    % draw roi
    bw = roipoly;
    ix = find(bw);
    if isempty(ix); break; end
    ind{end+1} = ix;
    
    % pull data
    hImage = imhandles(gca);
    cdata = get(hImage,'cdata');

    % do calculations
    [B,BINT,R2,RINT,STATS] = regress(cdata(ix),ones(size(ix)));
    n = numel(ix);
    m = B;
    ci95 = B-BINT(1);
    stats = [n m std(cdata(ix)) ci95];
    disp(['n / mean / std / ci95   :   ' num2str(stats)])
    
    % calm CPU
    pause(0.05);
    
end
set(h,'Position', siz);
if nargout<1; clear ind; end
