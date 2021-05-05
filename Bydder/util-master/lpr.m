function flag = lpr(filename)
%function lpr(filename)
% Exports current Figure to 'image.tif' on Desktop.
% If filename is supplied, writes to filename(.tif)
% or another recognized format (.png .jpg .bmp).

% no figures open!
if isempty(get(0,'currentfigure'))
    error('No figures currently open');
end

% default filename and image format
if nargin==0
	filename = 'image';
end
if isnumeric(filename)
    filename = num2str(filename);
end

% remove bad terminating characters
ok = 'abcdefghijklmnopqrstuvwxyz_0123456789';
while ~ismember(lower(filename(end)),ok)
    filename = filename(1:end-1);
    if isempty(filename); error('Filename not valid.'); end
end

% image file format
[pathstr,filestr,ext] = fileparts(filename);
switch ext
    case '';      format = 'tif'; ext = '.tif'; % default
    case '.png';  format = 'png';
    case '.jpg';  format = 'jpeg';
    case '.jpeg'; format = 'jpeg';
    case '.tif';  format = 'tif';
    case '.tiff'; format = 'tiff';
    case '.bmp';  format = 'bmp';
    otherwise; error('Image file format "%s" not recognized.',ext)
end

% default directory
if isempty(pathstr)
    pathstr = '~/Desktop';
end

% filename to write
filename = strcat(pathstr,filesep,filestr,ext);
if nargout
    flag = 0;
else
    fprintf('%s writing to %s\n',mfilename,filename);
end

% proceed
%Resolution = 1000; % default resolution (dots/inch?)
%hgexport(gcf, filename, hgexport('factorystyle'), 'Format', format, 'Resolution', 1000);
saveas(gcf,filename,format);
