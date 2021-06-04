classdef HamiltonLiver2011_FWspectrum < Custom_FWspectrum
    %Hamilton G, et al. NMR Biomed. 24(7):784-90, 2011. PMID: 21834002
    methods
        function this = HamiltonLiver2011_FWspectrum(varargin)
            this=this@Custom_FWspectrum(varargin{:});
            this.fatCS = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29];
            this.realAmps = [88, 642, 58, 62, 58, 6, 39, 10, 37] / 1000;
        end
    end
    
end

    
    
