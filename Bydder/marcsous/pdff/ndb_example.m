%% fit phantom data from the paper
%
% https://doi.org/10.1016/j.mri.2011.07.004
%

load PHANTOM_NDB_PAPER.mat

[nx ny ne] = size(data);
params = zeros(nx,ny,5);

% loop over x and y dimensions
for x = 1:nx
    
    ims(params,[],{'B0 (Hz)','R2* (1/s)','FF (%)','PHI (rad)','NDB'}); drawnow
    
    for y = 1:ny
        
        if any(data(x,y,:))
            params(x,y,:) = struct2array(ndbest(te,data(x,y,:),Tesla,H2O));
        end
        
    end
end
