classdef Berglund2012_FWspectrum < Custom_FWspectrum
    properties
    CL=17.4
    P2U=0.2
    UD=2.6
    end
    
    methods
        function this = Berglund2012_FWspectrum(varargin)
            this=this@Custom_FWspectrum(varargin{:});
            this.fatCS = [0.90, 1.30, 1.59, 2.03, 2.25, 2.77, 4.1, 4.3, 5.21, 5.31];
        end
        
        function amps=setFattyPeaks(this)
            CL = this.CL;
            P2U = this.P2U;
            UD = this.UD;
            amps = 0;
            switch this.UFV
                case 0
                    amps = [9, (CL - 4) * 6 + UD * (2 * P2U - 8), 6, 4 * (UD * (1 - P2U)), 6, 2 * UD * P2U, 2, 2, 1, 2 * UD];
                    amps = amps / sum(amps) ;
                
                case 1
                    % Unkwown Fat variable : UD
                    % F1 = 9A+6(CL-4)B+6C+6E+2G+2H+I
                    % F2 = (2P2U-8)B+4(1-P2U)D+2P2UF+2J
                    amps = [[9, 6 * (CL - 4), 6, 0, 6, 0, 2, 2, 1, 0];
                        [0, 2 * P2U - 8, 0, 4 * (1 - P2U), 0, 2 * P2U, 0, 0, 0, 2]];
                
                case 2
                    % Unkwown Fat variables : UD, P2U
                    % F1 = 9A+6(CL-4)B+6C+6E+2G+2H+I
                    % F2 = -8B+4D+2J
                    % F3 = 2B-4D+2F
                    amps = [[9, 6 * (CL - 4), 6, 0, 6, 0, 2, 2, 1, 0];
                        [0, -8, 0, 4, 0, 0, 0, 0, 0, 2];
                        [0, 2, 0, -4, 0, 2, 0, 0, 0, 0]];
                    
                case 3
                    % Unkwown Fat variables : UD, P2U, CL
                    % F1 = 9A-24B+6C+6E+2G+2H+I
                    % F2 = -8B+4D+2J
                    % F3 = 2B-4D+2F
                    % F4 = 6B
                    amps = [[9, -24, 6, 0, 6, 0, 2, 2, 1, 0];
                    [0, -8, 0, 4, 0, 0, 0, 0, 0, 2];
                    [0, 2, 0, -4, 0, 2, 0, 0, 0, 0];
                    [0, 6, 0, 0, 0, 0, 0, 0, 0, 0]];
                otherwise
                    error("Ukwnown Fat variable should be between 0 and 3")
            end
        
    end
    
    end
end
