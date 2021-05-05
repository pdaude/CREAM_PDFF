function [left, right, up, down, above, below] = getIndexImages(nx, ny, nz)
%GETINDEXIMAGES Gets index images
%   Gets index images

left = zeros(nz, ny, nx,'logical');
left(:, :, 1:end-1) = 1;
right = zeros(nz, ny, nx,'logical');
right(:, :, 2:end) = 1;
down = zeros(nz, ny, nx,'logical');
down(:, 1:end-1, :) = 1;
up = zeros(nz, ny, nx,'logical');
up(:, 2:end, :) = 1;
below = zeros(nz, ny, nx,'logical');
below(1:end-1, :, :) = 1;
above = zeros(nz, ny, nx,'logical');
above(2:end, :, :) = 1;

left = permute(left,[3,2,1]);
right = permute(right,[3,2,1]);
up = permute(up,[3,2,1]);
down = permute(down,[3,2,1]);
above = permute(above,[3,2,1]);
below = permute(below,[3,2,1]);

left = left(:);
right = right(:);
up = up(:);
down = down(:);
above = above(:);
below = below(:);

end