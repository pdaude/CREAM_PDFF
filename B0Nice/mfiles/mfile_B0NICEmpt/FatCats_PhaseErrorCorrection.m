function [algoParams] = FatCats_PhaseErrorCorrection(algoParams)
%--------------
%make sure the peak of water-pixels 
    % global phase correction  
    algoParams.BW_label = ones(algoParams.matrix_size(1:3)); % including all pixels
    [algoParams] = CheckB0(algoParams);
    %-------------------------
    % divide B0 map into regions using phase gradient
    [algoParams] = DefineRegion_usingPhaseGradient(algoParams);
    [algoParams] = CheckB0(algoParams);
    
    if algoParams.index_B0(1) == 1 && (algoParams.index_B0(3)-algoParams.index_B0(2)) > 2;
        % divide B0 map into regions using the difference of FF between
        % mag- and phase_based
        [algoParams] = DefineRegion_usingDiffFF(algoParams);
        [algoParams] = CheckB0(algoParams);
    end
%---------------------------------------------
