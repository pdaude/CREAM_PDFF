classdef Leporq2014SAT_FWspectrum < Hamiltonndb2011_FWspectrum
    %Leporq B, et al. NMR Biomed. 2014; 27: 1211–1221 DOI:10.1002/nbm.3175
    methods
        function this = Leporq2014SAT_FWspectrum(varargin)
            this = this@Hamiltonndb2011_FWspectrum(varargin{:});
            this.UFV= 0;
            this.NDB=2.53;
            this.NMIDB = 0.94;
            this.CL = 17.47;
            this.fatCS = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29];
            
        end
    end
    
end
