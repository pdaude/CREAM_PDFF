classdef ISMRMChallenge_FWspectrum < Custom_FWspectrum
    %Hamilton G, et al. NMR Biomed. 24(7):784-90, 2011. PMID: 21834002
    methods
        function this = ISMRMChallenge_FWspectrum(varargin)
            this=this@Custom_FWspectrum(varargin{:});
            this.fatCS = [0.90, 1.30, 2.1, 2.76, 4.31, 5.3];
            this.realAmps = [87, 693, 128, 4, 39, 48]/1000;
        end
    end
    
end

    
    
