function [species,FWSpectrumOutput] = setupModelParams(mPar,varargin) 
%setupModel Params is 
    opts.clockwisePrecession = false; 
    opts.temperature= []; 
    
    % varargin handling (must be option/value pairs)
    for k = 1:2:numel(varargin)
        if k==numel(varargin) || ~ischar(varargin{k})
        error('''varargin'' must be option/value pairs.');
        end
        if isfield(opts,varargin{k})
           opts.(varargin{k}) = varargin{k+1};
        end
    end
    
    defaultmPar = defaultModelParams();
    defaultfields=fieldnames(defaultmPar);
    
    for n= 1:numel(defaultfields)
        strfield=string(defaultfields(n));
         if ~isfield(mPar ,defaultfields(n))
             mPar.(strfield)=defaultmPar.(strfield);
         end
    end
    
    if ~isempty(opts.temperature) %Temperature dependence according to Hernando 2014
         mPar.watCS = 1.3 + 3.748 - .01085 * opts.temperature;  % Temp in [°C]
    end
    
     %List of all FWSpectrumClasses
     [dirpath,name,ext] = fileparts(mfilename('fullpath'));
     list_dirFWSpectrum=dir(fullfile(dirpath,'FWSpectrumClasses','@*_FWspectrum'));
     FWSpectrumClasses=cellfun(@(x) x(2:end),{list_dirFWSpectrum.name},'UniformOutput',false);
     FWSpectrumClassesname= cellfun(@(x) x(2:strfind(x,'_FWspectrum')-1),{list_dirFWSpectrum.name},'UniformOutput',false);
     FWSpectrumClass=FWSpectrumClasses(strcmp(FWSpectrumClassesname,mPar.model));
     if isempty(FWSpectrumClass)
        error('Unknown spectrum model :  \n Available models are %s',string(join(FWSpectrumClassesname,', ')))
     end
     X = fieldnames(mPar);
     Y = struct2cell(mPar);
     C = [X,Y].';
     FWSpectrumOutput=feval(string(FWSpectrumClass),C{:});
     species=FWSpectrumOutput.ISMRMformat();
    
end

function defaultsmPar=defaultModelParams()
defaultsmPar.model = 'Berglund2012';
defaultsmPar.watCS = 4.7;
defaultsmPar.gyro = 42.577478;
end

