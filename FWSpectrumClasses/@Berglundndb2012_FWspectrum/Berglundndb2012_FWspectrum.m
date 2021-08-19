classdef Berglundndb2012_FWspectrum < Custom_FWspectrum
    properties
    NDB=0;
    NMIDB=0;
    CL=0;
    end
    
    methods
        function this = Berglundndb2012_FWspectrum(varargin)
            this=this@Custom_FWspectrum(varargin{:});
            this.fatCS = [0.90, 1.30, 1.59, 2.03, 2.25, 2.77, 4.1, 4.3, 5.21, 5.31];
        end
        
        function amps=setFattyPeaks(this)
            CL = this.CL;
            ndb = this.NDB;
            nmidb = this.NMIDB;
            amps = 0;
            switch this.UFV
                case 0
                    amps = [9, ((CL - 4) * 6) - (ndb * 8) + nmidb * 2, 6, (ndb - nmidb) * 4, 6, nmidb * 2, 2, 2, 1, 2 * ndb];
                    amps = amps / sum(amps);
                
                case 1
                    %Unkwown Fat variable : ndb
                    % F1 = 9A+(6(CL-4)+2nmidb)B+6C-4(nmidb)D+6E+2(nmidb)F+4(G+H)+I
                    % F2 = -8B+4D+2J
                    amps = [[9, 6 * (CL - 4)+2*nmidb, 6, nmidb*-4, 6, nmidb * 2, 2, 2, 1, 0];
                             [0,- 8, 0,4, 0, 0, 0, 0, 0, 2]];
                
                case 2
                    % Unkwown Fat variable : ndb, nmidb
                    % F1 = 9A+(6(CL-4)B+6C+6E+4(G+H)+I
                    % F2 = -8B+4D+2J
                    % F3 = 2B -4D+2F
                    amps = [[9, 6 * (CL - 4), 6,0, 6,0, 2, 2, 1, 0];
                             [0, - 8, 0, 4, 0, 0, 0, 0, 0, 2];
                             [0, 2, 0, -4, 0, 2, 0, 0, 0, 0]];
                    
                case 3
                    % Unkwown Fat variable = ndb, nmidb and CL
                    % F1 = 9A-24B+6C+6E+4(G+H)+I
                    % F2 = -8B+4D+2J
                    % F3 = 2B -4D+2F
                    % F4 = 6B
                    amps = [[9, -24, 6, 0, 6, 0, 2, 2, 1, 0];
                             [0, - 8, 0, 4, 0, 0, 0, 0, 0, 2];
                             [0, 2, 0, -4, 0, 2, 0, 0, 0, 0];
                             [0, 6, 0, 0, 0, 0, 0, 0, 0, 0]];
                otherwise
                    error("Ukwnown Fat variable should be between 0 and 3")
            end
        
    end
    
    end
end
