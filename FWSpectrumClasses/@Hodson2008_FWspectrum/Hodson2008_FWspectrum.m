classdef Hodson2008_FWspectrum < Berglundndb2012_FWspectrum
    %Hodson L, et al. Prog Lipid Res. 2008;47:348-380. PMID: 18435934 DOI: 10.1016/j.plipres.2008.03.003
    methods
        function this = Hodson2008_FWspectrum(varargin)
            this = this@Berglundndb2012_FWspectrum(varargin{:});
            this.UFV= 0;
            this.NDB=2.69;
            this.NMIDB = 0.58;
            this.CL = 17.29;
            this.fatCS = [0.90, 1.30, 1.59, 2.03, 2.25, 2.77, 4.1, 4.3, 5.21, 5.31];
            
        end
    end
    
end
