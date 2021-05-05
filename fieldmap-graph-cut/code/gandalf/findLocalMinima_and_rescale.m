function [nMinimaPerVoxel, costLocalMinimaRescale, indexLocalMinima] = findLocalMinima_and_rescale(residual, masksignal, options)

    nMaxMinimizers = options.nMaxMinimizers;
    [Numlayers, nVoxel_Y, nVoxel_X] = size(residual);
    residual_islocalminima = zeros(Numlayers, nVoxel_Y, nVoxel_X);
    
    %% check for local minimizers
    for Y = 1:nVoxel_Y
        for X = 1:nVoxel_X
            if masksignal(Y, X) ~= 0
                residual_islocalminima(:, Y, X) = islocalmin(residual(:, Y, X), 'MinSeparation', options.minDistance);
            end
        end
    end
    
    nMinimaPerVoxel = squeeze( sum( residual_islocalminima, 1));
    nMinimaPerVoxel = nMinimaPerVoxel .* masksignal;
%     nMaxNodesPerVoxel = max(nMinimaPerVoxel(:));
    
    costLocalMinima = zeros(nVoxel_Y, nVoxel_X, nMaxMinimizers);
    indexLocalMinima = zeros(nVoxel_Y, nVoxel_X, nMaxMinimizers);
    
    %% Extract all sets of local minimizers
    for Y = 1:nVoxel_Y
        for X = 1:nVoxel_X
            rem = 1;
            for Z = 1:Numlayers
                if residual_islocalminima(Z, Y, X) == 1
                    costLocalMinima(Y, X, rem) =  residual(Z, Y, X);                  
                    indexLocalMinima(Y, X, rem) = Z;
                    rem = rem + 1;
                end
            end
        end
    end
    
    %% rescale all sets of local minimizers to $intervalVector
    costLocalMinimaRescale = zeros(nVoxel_Y, nVoxel_X, nMaxMinimizers);
    for Y = 1:nVoxel_Y
        for X = 1:nVoxel_X
            if masksignal(Y, X) ~= 0 && (nMinimaPerVoxel(Y, X) > 1)
                tmp_diff = diff(costLocalMinima(Y, X, 1:nMinimaPerVoxel(Y, X)));
                tmp_diff = sum(tmp_diff);
                if tmp_diff ~= 0
                    costLocalMinimaRescale(Y, X, 1:nMinimaPerVoxel(Y, X)) = scale_array2interval(costLocalMinima(Y, X, 1:nMinimaPerVoxel(Y, X)), options.rescale);
                else
                    costLocalMinimaRescale(Y, X, 1:nMinimaPerVoxel(Y, X)) = 1;
                end
            end
        end
    end
    
end