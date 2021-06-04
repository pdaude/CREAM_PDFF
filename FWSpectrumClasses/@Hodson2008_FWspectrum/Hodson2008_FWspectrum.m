classdef Hodson2008_FWspectrum < Hamiltonndb2011_FWspectrum
    %Hodson L, et al. Prog Lipid Res. 2008;47:348-380. PMID: 18435934 DOI: 10.1016/j.plipres.2008.03.003
    methods
        function this = Hodson2008_FWspectrum(varargin)
            this = this@Hamiltonndb2011_FWspectrum(varargin{:});
            this.UFV= 0;
            this.NDB=2.69;
            this.NMIDB = 0.58;
            this.CL = 17.29;
            this.fatCS = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29];
            
        end
    end
    
end
