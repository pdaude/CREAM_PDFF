% Function: simuGenericTE
%
% Description: Generate echo times with uniform spacing meeting hardware
% constraints (TEmin and dTemin) and echo type (IN/OPP; IDEAL;MINIMUM)
% 
% Parameters:
% Input: 
% Mandatory inputs
%  - NTE: number of echo times 
%   - B0 : field strength in Tesla
%  - type : Type of echo (IN/OPP; IDEAL)
% Optionnal inputs
%   - dTEmin : Minimum echo spacing constraint by hardware 
%   - TEmin : Minimum echo time constraint by hardware 
%   - NTEmax : number of 2 pi echo times to assess constraints 
%
% Returns: 
% - TE : echo times array with uniform spacing meeting hardware constraints and echo type

% Author: Pierre Daudé
function TE = simuGenericTE(NTE,B0,type,varargin)

opts.dTEmin=1.37e-3; % VIDA WIP Lausanne sequence
opts.TEmin=1.09e-3; %VIDA WIP Lausanne sequence
opts.NTEmax = 20; 

CSFW=4.7-1.3; % ppm
gyro=42.57747892;
larmor=B0*gyro; %(MHz)

for k = 1:2:numel(varargin)
    if k==numel(varargin) || ~ischar(varargin{k})
        error('''varargin'' must be option/value pairs.');
    end
    if ~isfield(opts,varargin{k})
        error('''%s'' is not a valid option.',varargin{k});
    end

    opts.(varargin{k}) = varargin{k+1};
end

% Assess Multi Echo or not 
if NTE < 2
error('''NTE'' is less than to 2');
end


% Assess echo type
switch type
    case 'IN/OPP'
        
        % angle=pi+kpi k in [0;NTE-1]
        % angle= angle_min + kpi k in [0;NTE-1]
    
        % Calculate INOPP dTE which corresponds to theta=pi+k2pi
        Ntype= 2;
        angle2pik_n=pi;
        
    case 'IDEAL'
        
        Ntype= NTE;
        %angle2pik_n=(0.5+(2*(0:Ntype-1)/Ntype))*pi;%step inside 2*pi : k/N*2pi k in [1:Ntype-1]  
    
    case 'MINIMUM'
  
        TE=opts.TEmin+(0:1:NTE-1)*opts.dTEmin;
        return
        
    otherwise
        error('type should be IDEAL and IN/OPP or MINIMUM');
end


% Determine the minimum delta TE which repects the sequence type (IN/OPP or
% IDEAL and the Hardware (opts.dTEmin)

p=primes(floor(Ntype/2)); % primes number without 1
k=1:Ntype-1;
for i=p
   if mod(Ntype,i)==0  % SI i factorielle de Ntype 
       k=k(mod(k,i)~=0) ; %Remove all  k multiple de i 
   end
end

step2pik= 2*(0:opts.NTEmax-1)*pi; %2*pi step  for IN/OPP TEs
step2pik_n= 2*k/Ntype*pi;  %step inside 2*pi : k/N*2pi k in [1:Ntype-1]  

[step2piK,step2piK_N] = meshgrid(step2pik,step2pik_n);

dTEs=reshape((step2piK+step2piK_N)/(2*pi*larmor*CSFW),[numel(step2piK),1]);

idx_min_dTE=find(dTEs-opts.dTEmin>=0);

if numel(idx_min_dTE)==0
        error('''dTEmin'' is too big. Please specified a NTEmax bigger');
end

if idx_min_dTE(1)>= Ntype
        warning(' With this ''dTEmin'', your step  is superior to 2pi');
end

step_min=dTEs(idx_min_dTE(1));

% Determine the minimum TE which repects the sequence type (IN/OPP or
% IDEAL and the Hardware (opts.TEmin)
switch type
    case 'IN/OPP'
        
        [step2piK,angle2piK_N] = meshgrid(step2pik,angle2pik_n);

        TEs=reshape((step2piK+angle2piK_N)/(2*pi*larmor*CSFW),[numel(step2piK),1]);

        idx_min_TE=find(TEs-opts.TEmin>=0);

        if numel(idx_min_TE)==0
                error('''TEmin'' is too big. Please specified a NTEmax bigger');
        end

        if idx_min_TE(1)>= Ntype
                warning(' With this ''dTEmin'', your TEmin is superior to 2pi');
        end

        TE_min=TEs(idx_min_TE(1));
        
    case 'IDEAL'
        TE_min=opts.TEmin;
        
    otherwise
        error('type should be IDEAL and IN/OPP');
end

TE = TE_min+(0:NTE-1)*step_min;

end


