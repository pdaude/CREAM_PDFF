classdef Bydder2011_FWspectrum < Hamiltonndb2011_FWspectrum
    %Hamilton G, et al. NMR Biomed. 24(7):784-90, 2011. PMID: 21834002
    methods
        function this = Bydder2011_FWspectrum(varargin)
            this = this@Hamiltonndb2011_FWspectrum(varargin{:});
            if this.NDB==0 
                this.NDB=2.5;
            end
            this.fatCS = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29];
            this.NMIDB = 0.093*this.NDB.^2;
            this.CL = 16.8 + 0.25*this.NDB;
            
        end
    end
    
end

    
    