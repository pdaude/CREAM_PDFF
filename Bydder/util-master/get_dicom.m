function [data,TE,SL,FA,TR,TI,BV,LF,SN,AN,dcm,pathname] = get_dicom(token,pathname,maxno)
%
% Usage: [data,TE,SL,FA,TR,TI,BV,LF,SN,AN,dcm,pathname] = get_dicom(token,pathname,maxno);
%
% Locates all DICOM files in the path specified by pathname and filtered by token.
% Recurses into all subdirectories.
%
% Inputs:
%  token is an identifier to parse the DICOM header:
%   -if numeric then we assume it's a seriesNumber
%   -if string then we assume it's SeriesDescription
%   -if struct then we assume it's returned from loaddcmdir.m
%  pathname is the path to search in (searches recursively)
%  maxno is the max. number of files to read in
%
% Outputs:
%  data is a high dimensional array sorted by parameters

% parallel option (only works in later matlab versions)
try; gcp; end

%% handle arguments
if ~exist('token','var')
    token = []; % empty string also acceptable
end
if ~exist('pathname','var') | isempty(pathname)
    pathname = uigetdir();
    if isequal(pathname,0); return; end
end
if isa(token,'struct')
    pathname = token.Path;
end

% convert file separators for local machine
pathname(pathname=='\' | pathname=='/') = filesep;
if pathname(end) ~= filesep
    pathname(end+1) = filesep;
end
if ~exist(pathname,'dir')
    error(['pathname does not exist: ' pathname])
end
if ~exist('maxno','var') | isempty(maxno)
    maxno = Inf;
end
disp([mfilename '(): ' pathname])

% parse token argument
if isa(token,'struct')
    disp([mfilename '(): reading in files specified by loaddcmdir(). Found 0   '])
    info = cell2struct(token.Images,'name',numel(token.Images));
else
    if isempty(token)
        disp([mfilename '(): Reading in all files. Found 0   '])
    elseif isa(token,'char')
        disp([mfilename '(): Parsing SeriesDescription with ''' token '''. Found 0   '])
    elseif isa(token,'numeric')
        disp([mfilename '(): Parsing SeriesNumber with ''' num2str(token) '''. Found 0   '])
    end
    info = get_files(pathname); % recurse into directories
end

% allow read-in from multiple series (false if filtering by series no.)
if isa(token,'numeric') & ~isempty(token)
    ignore_SN = false;
else
    ignore_SN = true;
end

%% read in data

% create file counter using temp file (workaround for parfor)
tmpfile = [tempname '.delete.me'];
fid = fopen(tmpfile,'w');
fclose(fid);

parfor j = 1:numel(info) % if there is an error here, change parfor to for

    % current file
    filename = [pathname info(j).name];
    
    % test for maxno and whether it's a DICOM file
    if update_counter(tmpfile,maxno) & isdicom(filename)
 
        % read header
        head = dicominfo(filename);
        
        % filter by property
        if isempty(token) | isa(token,'struct')
            found = true;
        elseif isa(token,'char')
            found = ~isempty(strfind(head.SeriesDescription,token));
            %keyboard
            %found = ~isempty(strfind(head.PatientName,token));
        elseif isa(token,'numeric')
            found = head.SeriesNumber==token;
        else
            found = false;
        end

        % store properties
        if found

            % update count
            update_counter(tmpfile);
            
            % read straightforward tags
            TE(j) = head.EchoTime;
            SL(j) = head.SliceLocation;
            FA(j) = head.FlipAngle;
            TR(j) = head.RepetitionTime;
            IN(j) = head.InstanceNumber;
            AN(j) = head.AcquisitionNumber;
            SN(j) = head.SeriesNumber;
            LF(j) = head.ImagingFrequency;
            
            % real/imag flags: 0=mag 1=phase 2=real 3=imag 
            RI(j) = 0; % default if no tag is present
            if isfield(head,'Private_0043_102f') % GE
                RI(j) = head.Private_0043_102f(1);
            end
            if isfield(head,'Private_0051_1016') % Siemens
                % NEEDS CHECKING FOR R/I - MAY BE I/Q ?
                RI(j) = strfind('MPRI',upper(head.Private_0051_1016(1)))-1;
            end

            % inversion tag not always present
            if isfield(head,'InversionTime')
                TI(j) = head.InversionTime;
            end
            
            % b-values not stored in standard tags... not working for DTI
            if isfield(head,'Private_0043_1039')
                %tensor = [head.Private_0019_10bb ...
                %          head.Private_0019_10bc ...
                %          head.Private_0019_10bd];
                
                % weird GE thing to add 1e9 except when b=0
                BV(j) = max(0,head.Private_0043_1039(1)-1e9);              
            end
            
            % store images in cell array (less memory)
            data{j} = dicomread(filename);
            
            % store header
            dcm{j} = head;

        end % is found
        
    end % is dicom
    
end

% clean up
delete(tmpfile);

% check for empty
if ~exist('data')
    data = [];
    TE = [];
    SL = [];
    FA = [];
    TR = [];
    TI = [];
    BV = [];
    LF = [];
    SN = [];
    AN = [];
    dcm = [];
    return;
end

% count non-empty data cells
index = find(~cellfun('isempty',data));
count = numel(index);

% parfor counter increments in multiples of ncpu, trim excess
if isfinite(maxno)
    count = maxno;
    index = index(1:maxno);
end

% parfor counter is flaky so re-display actual count
disp(sprintf('\b\b\b\b\b%-4d',count));

% convert cell arrays to matrices
[nx ny] = size(data{index(1)});
for j = 1:count
    temp = size(data{index(j)});
    if any(temp - [nx ny])
        error('size mis-match %i (expecting %ix%i found %ix%i)',j,nx,ny,temp(1),temp(2));
    end
end
dcm = dcm(index);
data = data(index);
data = cell2mat(data);

data = reshape(data,nx,ny,count); % image data
TE = TE(index); % echo time
SL = SL(index); % slice location
FA = FA(index); % flip angle
TR = TR(index); % repetition time
IN = IN(index); % instance number
SN = SN(index); % series number
RI = RI(index); % rawdata type (real/imag)
LF = LF(index); % larmor freq
AN = AN(index); % acquisition no. (repetitions)
if exist('TI','var'); TI = TI(index); else TI = []; end % inversion time
if exist('BV','var'); BV = BV(index); else BV = []; end % bvalue 

% unique properties
[uTE,~,kte] = unique(TE);
[uSL,~,ksl] = unique(round(SL*1e3)/1e3); % round floats (sometimes values have jitter)
[uFA,~,kfa] = unique(FA);
[uTR,~,ktr] = unique(TR);
[uSN,~,ksn] = unique(SN);
[uTI,~,kti] = unique(TI);
[uBV,~,kbv] = unique(BV);
[uRI,~,kri] = unique(RI);
[uLF,~,klr] = unique(round(LF*1e3)/1e3); % round floats (sometimes values have jitter)
[uAN,~,kan] = unique(AN);

nTE = numel(uTE);
nSL = numel(uSL);
nFA = numel(uFA);
nTR = numel(uTR);
nSN = numel(uSN);
nTI = numel(uTI);
nBV = numel(uBV);
nRI = numel(uRI);
nLF = numel(uLF);
nAN = numel(uAN);
if nTI==0; nTI = 1;	kti = ones(count,1); end
if nBV==0; nBV = 1;	kbv = ones(count,1); end
if nLF~=1; error('More than one field strength present!'); end
if nAN>1 && nSN==1 % this handles the case where we vary parameters over the repetition loop
    disp([mfilename '(): Detected ' num2str(nAN) ' repetitions in 1 series. Combining into 1 repetition.'])    
    nAN = 1;
    kan = ones(size(kan));
end
if ignore_SN && nSN>1
    disp([mfilename '(): Combining ' num2str(nSN) ' different series. Specify SeriesNumber to avoid this.'])
    for j = 1:nSN
        disp(['  ' num2str(uSN(j)) '  ' dcm{j}.SeriesDescription])
    end
    nSN = 1;
    uSN = 1;
    ksn = ones(size(ksn));
end

% catch size errors
expected_count = nTE*nSL*nFA*nTR*nTI*nBV*nRI*nSN*nAN;

if count~=expected_count
    warning('%s() wrong number of images (found %i but expecting %i)',mfilename,count,expected_count)
    disp(['CHECK THESE!! echos=' num2str(nTE) ...
        ' slices=' num2str(nSL) ...
        ' flips=' num2str(nFA) ...
        ' TRs=' num2str(nTR) ...
        ' TIs=' num2str(nTI) ...
        ' Bs=' num2str(nBV) ...
        ' series=' num2str(nSN) ' (ignore=' num2str(ignore_SN) ')'...
        ' real/imag=' num2str(nRI) ... 
        ' aquisitions=' num2str(nAN)])
    disp('Type "return" to ignore (expect errors!) or dbquit to cancel.')
    keyboard
end

% sort by property
temp = zeros(nx,ny,nSL,nTE,nFA,nTR,nTI,nBV,nSN,nRI,nAN,'single');
dcmtemp = cell(nSL,nTE,nFA,nTR,nTI,nBV,nSN,nRI,nAN);
mask = zeros(size(dcmtemp)); % for debugging - counts the no. of images in each slot
for j = 1:count
    temp(:,:,ksl(j),kte(j),kfa(j),ktr(j),kti(j),kbv(j),ksn(j),kri(j),kan(j)) = data(:,:,j);
    dcmtemp(ksl(j),kte(j),kfa(j),ktr(j),kti(j),kbv(j),ksn(j),kri(j),kan(j)) = dcm(j);
    mask(ksl(j),kte(j),kfa(j),ktr(j),kti(j),kbv(j),ksn(j),kri(j),kan(j)) = ...
        mask(ksl(j),kte(j),kfa(j),ktr(j),kti(j),kbv(j),ksn(j),kri(j),kan(j))+1;
end
data = temp;
dcm = dcmtemp;
TE = uTE;
SL = uSL;
FA = uFA;
TE = uTE;
TR = uTR;
SN = uSN;
BV = uBV;
LF = uLF;
AN = uAN;
if any(mask~=1)
    disp('Something is wrong... duplicates or missing data - need to check')
    keyboard
end
clear temp dcmtemp mask;

%% handle real/imag/mag/phase (not perfect)

% scale phase from int into float range (0 to 2pi)
if isequal(dcm{1}.Manufacturer,'SIEMENS')
    phase_scale = single(2*pi/4095);
elseif isequal(dcm{1}.Manufacturer,'GE')
    phase_scale = single(2*pi/1000);
else
    phase_scale = 1;
    if any(uRI==1)
        warning('get-dicom() phase_scale not defined for phase images')
    end
end

if nRI==1
    % 1=phase
    if isequal(uRI,1); data = data*phase_scale-pi; end
elseif nRI==2 & isequal(uRI,[0 1])
    % 0=mag and 1=phase
    data = data(:,:,:,:,:,:,:,:,:,1).*exp(i*(data(:,:,:,:,:,:,:,:,:,2)*phase_scale-pi));
elseif nRI==2 & isequal(uRI,[2 3])
    % 2=real and 3=imag
    data = complex(data(:,:,:,:,:,:,:,:,:,1),data(:,:,:,:,:,:,:,:,:,2));
elseif nRI==4 & isequal(uRI,[0 1 2 3])
    % ignore mag/phase (less reliable than real/imag)
    data = complex(data(:,:,:,:,:,:,:,:,:,3),data(:,:,:,:,:,:,:,:,:,4));
    disp([mfilename '(): real/imag/mag/phase present - returning real/imag'])
else
    disp([mfilename '(): not sure about real/imag/phase/mag - need to handle this'])
    keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recursively fetch file info from a directory %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function info = get_files(pathname)

% contents of starting directory
info = dir(pathname);

j = 1;
while j <= numel(info)
    % recurse into subdirectories
    if info(j).isdir
        % skip subdirectories with a leading '.'
        if info(j).name(1) ~= '.'
            temp = dir([pathname info(j).name]);
            % prepend path (except for '.')
            for k = 1:numel(temp)
                if temp(k).name(1) ~= '.'
                    temp(k).name = [info(j).name filesep temp(k).name];
                end
            end
            % append contents of subdirectory
            info = [info;temp];
        end
        % delete directory from list
        info(j) = [];
    else
        % skip past files
        j = j + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect if file is DICOM (borrowed from old MATLAB version) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = isdicom(filename)
%ISDICOM    Determine if a file is probably a DICOM file.
%    TF = ISDICOM(FILENAME) returns true if the file in FILENAME is
%    probably a DICOM file and FALSE if it is not.

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/06/15 20:10:36 $

% MB exclude DICOMDIR, which otherwise returns true
if strfind(filename,'DICOMDIR')
    tf = false;
    return;
end

% Open the file.
fid = fopen(filename, 'r');
if (fid < 0)
    %warning('Image:isdicom:fileOpen', ...
    %    'Couldn''t open file for reading: %s',filename)
    tf = false;
    return % MB just return, don't complain
end

% Get the possible DICOM header and inspect it for DICOM-like data.
header = fread(fid, 132, 'uint8=>uint8');
fclose(fid);

if numel(header)<132 % MB exclude easy pickings

    % It's too small
    tf = false;
    
elseif isequal(char(header(129:132))', 'DICM')
    
    % It's a proper DICOM file.
    tf = true;

else
    % MB heuristic is not good enough - exclude all.
    %
    % Use a heuristic approach, examining the first "attribute".  A
    % valid attribute will likely start with 0x0002 or 0x0008.
    %group = typecast(header(1:2), 'uint16');
    %if (isequal(group, uint16(2)) || isequal(swapbytes(group), uint16(2)) || ...
    %        isequal(group, uint16(8)) || isequal(swapbytes(group), uint16(8)))
    % 
    %    tf = true;
    %else
        tf = false;
    %end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         tricky counter for parfor loop           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flag = update_counter(tmpfile,maxno)
% Warning: "clever"
% 1 input = increment counter
% 2 inputs = display counter and test for maxno
if nargin==1 

    % only proceed if file is not in use by another lab
    fid = fopen(tmpfile,'a');
    while fid==-1
        pause(rand*1e-2);
        fid = fopen(tmpfile,'a');
    end
 
    % append another 1 to counter file
    fwrite(fid,1);
    fclose(fid);
    
else
    
    % tmpfile is a vector (1 entry per image)
    fid = fopen(tmpfile,'r');
    count = numel(fread(fid));
    fclose(fid);
    
    % counter status
    flag = count<maxno;
    
    % display counter
    if flag
        % always use 4 characters
        disp(sprintf('\b\b\b\b\b%-4d',count));
    end

end
