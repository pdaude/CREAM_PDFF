function map = rainbow()
%RAINBOW    rainbow color map
%   RAINBOW() returns an 256-by-3 matrix containing a "rainbow" colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(rainbow)
%
%   See also HOT, HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, 
%   COLORMAP, RGBPLOT.

%   E. Brian Welch, 2014-05-21.

% parabola vertex form with vertex at (h,k)
% y = a*(x-h)^2 + k
h = 29;
k = 0.35;

% also passes through (1,0)
a = -k/(1-h)^2;
x = [1:57]';
R_parabola = a*(x-h).^2 + k;

R = [ R_parabola; zeros(146-57,1); linspace(0,1,192-146)' ; ones(256-192,1) ];
G = [ zeros(57,1); linspace(0,1,102-57)' ; ones(192-102,1) ; linspace(1,0,256-192)' ];
B = [ linspace(0,1,53)' ; ones(102-53,1) ; linspace(1,0,147-102)'; zeros(256-147,1)];

map = [R G B];