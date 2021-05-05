function make_a_movie(data,CLIM,filename)

% Takes a 3D data set and writes an avi movie to 'movie.avi'
% -data is the 3D data set
% -CLIM=[CLOW CHIGH] is passed to imagesc (optional)
% -filename will write the file to filename.avi

if ~exist('filename','var')
    filename = 'movie';
end

% new way
M = VideoWriter([filename '.avi']);
open(M);

for j = 1:size(data,3)
    if nargin<3 || isempty(CLIM)
        imagesc(data(:,:,j));
    else
        imagesc(data(:,:,j),CLIM);
    end
    axis off;
    title(filename)
    drawnow
    writeVideo(M,getframe(gcf));
end
close(M);

disp(['Written file: ' pwd filesep filename '.avi'])