function H = getKernel(dx, dy, dz, sigma, use3D)
%GETKERNEL Gets kernel
%   Gets kernel
% HSIZE_x = 2*ceil(2*sigma/dx)+1;
% HSIZE_y = 2*ceil(2*sigma/dy)+1;
% HSIZE_z = 2*ceil(2*sigma/dz)+1;

if use3D
    
    [xx,yy,zz] = meshgrid(-ceil(2*sigma/dx):ceil(2*sigma/dx),-ceil(2*sigma/dy):ceil(2*sigma/dy),-ceil(2*sigma/dz):ceil(2*sigma/dz));
    
    dist = sqrt((dx*xx).^2 + (dy*yy).^2 + (dz*zz).^2);
    H = normpdf(dist,0,sigma);
    H = H./sum(H(:));
    
else
    
    [xx,yy] = meshgrid(-ceil(2*sigma/dx):ceil(2*sigma/dx),-ceil(2*sigma/dy):ceil(2*sigma/dy));
    
    dist = sqrt((dx*xx).^2 + (dy*yy).^2);
    H = normpdf(dist,0,sigma);
    H = H./sum(H(:));
    
end

end