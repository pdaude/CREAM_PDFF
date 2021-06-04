classdef Custom_FWspectrum
    properties
        UFV =0 ;
        fatCS = [];
        realAmps = [];
        watCS = 4.7;
        gyro = 42.577478;
    end
    methods
        function this =Custom_FWspectrum(varargin)
            for k = 1:2:numel(varargin)
                if k==numel(varargin) || ~ischar(varargin{k})
                    error('''varargin'' must be option/value pairs.');
                end
                if isprop(this,varargin{k})
                    this.(varargin{k}) = varargin{k+1};
                end
            end
        end
        function alpha = getAlpha(this)
            alpha=zeros(2+this.UFV,length(this.fatCS)+1,'double');
            alpha(1,1) = 1; % Water Component
            alpha(2:end,2:end) = this.setFattyPeaks();
            
        end
        
        function amps=setFattyPeaks(this)
            amps=this.realAmps;
        end
        
        function species = ISMRMformat(this)
            species(1).name = 'water';
            species(1).frequency = this.watCS;
            species(1).relAmps = 1;
            species(2).name = 'fat';
            species(2).frequency = this.fatCS;
            species(2).relAmps = this.setFattyPeaks();
        end
        
        
    end
end
