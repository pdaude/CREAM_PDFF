%% Function name: GANDALF
%%
%% Description: Fat-water separation using a single-step graphcut with unequally sampled space and penalized fieldmap jumps.
%%              + Variable Layer Graph Construction
%%              The whole 3D Volume is minimized in a single step.
%%
%% AUTHOR:      Christof Boehm <christof.boehm@tum.de>
%% AFFILIATION: Body Magnetic Resonance Research Group
%%              Department of Diagnostic and Interventional Radiology
%%              Technical University of Munich
%%              Klinikum rechts der Isar
%%              Ismaninger Str. 22, 81675 Muenchen
%% URL:         http://www.bmrrgroup.de
%%
%% Date created: May 14, 2018
%% Date last modified: September 7, 2018

function outParams = GANDALF(imDataParams, algoParams, VARPROparams)  
    
    % algoParams.range_r2star = VARPROparams.range_r2star;
    % algoParams.NUM_R2STARS = VARPROparams.NUM_R2STARS;
    INFTY = 100000000;
    
    %% Define and calculate key numbers of the input image
    [nVoxel_Y, nVoxel_X, nVoxel_Z, ~, ~] = size(imDataParams.images);

    costLocalMinimaRescale = VARPROparams.costLocalMinimaRescale;
    indexLocalMinimaRescale = VARPROparams.indexLocalMinimaRescale;
    masksignal = VARPROparams.masksignal;
    nMinimaPerVoxel = VARPROparams.nMinimaPerVoxel;
    
    nMaxNodesPerVoxel = max(nMinimaPerVoxel(:));
    nNodes = sum(nMinimaPerVoxel(:));
    nMeanfullVoxel = sum(masksignal(:));
    nVoxel = nVoxel_X * nVoxel_Y * nVoxel_Z;
    S_Node_ID = nNodes + 1;
    T_Node_ID = nNodes + 2;
    
    %% Calculate Noise Weighting
    MIP = get_echoMIP(imDataParams.images);
    NoiseWeighting = scale_array2interval(MIP, [0, 1]) .* masksignal;

    %% Create Node Indice Array
    t1 = tic;
    NodeIndicies = zeros(nVoxel_Y, nVoxel_X, nVoxel_Z, nMaxNodesPerVoxel);
    rem = 1;
    for Z = 1:nVoxel_Z
        for Y = 1:nVoxel_Y
            for X = 1:nVoxel_X
                for N = 1:nMinimaPerVoxel(Y, X, Z)
                   NodeIndicies(Y, X, Z, N) = rem;
                   rem = rem + 1;
                end
            end
        end
    end
    t1 = toc(t1);
    fprintf('Created Node Indice Array in %.2fs\n', t1);


    %% Create Intra-Column Edgelist
    t1 = tic;
    IntraColumnEdges = zeros((nNodes + nVoxel), 4);
    rem = 1;

    for Z = 1:nVoxel_Z
        for Y = 1:nVoxel_Y
            for X = 1:nVoxel_X
                if masksignal(Y, X, Z) ~= 0
                    for N = 1:nMinimaPerVoxel(Y, X, Z) - 1

                        node = NodeIndicies(Y, X, Z, N);
                        IntraColumnEdges(rem,1) = node;
                        IntraColumnEdges(rem,2) = node + 1;
                        IntraColumnEdges(rem,3) = costLocalMinimaRescale(Y, X, Z, N) * NoiseWeighting(Y, X, Z);
                        IntraColumnEdges(rem,4) = INFTY;
                        rem = rem + 1;

                        if N == 1                
                            IntraColumnEdges(rem, 1) = S_Node_ID;
                            IntraColumnEdges(rem, 2) = node;
                            IntraColumnEdges(rem, 3) = INFTY;
                            IntraColumnEdges(rem, 4) = 0;
                            rem = rem + 1;    
                        end

                        if N == (nMinimaPerVoxel(Y, X, Z) - 1)
                            IntraColumnEdges(rem, 1) = node+1;
                            IntraColumnEdges(rem, 2) = T_Node_ID;
                            IntraColumnEdges(rem, 3) = costLocalMinimaRescale(Y, X, Z, N + 1) * NoiseWeighting(Y, X, Z);
                            IntraColumnEdges(rem, 4) = 0;
                            rem = rem + 1;
                        end
                    end
                end
            end
        end
    end
    IntraColumnEdges = IntraColumnEdges(1:rem-1, :);
    t1 = toc(t1);
    fprintf('Created Intra-Column Edges in %.2fs\n', t1);

    %% Create Inter-Column Edgelist
    % Edges in X-Direction

    fprintf('Calculating Edges in X-Direction... ');

    t1 = tic;
    
    meanNodesPerVoxel = nNodes / nMeanfullVoxel;
    if nVoxel_Z == 1
        InterColumnEdges = zeros( ceil(3 * nMeanfullVoxel * meanNodesPerVoxel^2), 3); 
    else
        InterColumnEdges = zeros( ceil(2 * nMeanfullVoxel * meanNodesPerVoxel^3), 3); 
    end

    
    rem = 1;
    reverseStr = '';

    for Z = 1:nVoxel_Z
        for Y = 1:nVoxel_Y

            Prozent = ((((Z - 1) * nVoxel_Y) + Y) / (nVoxel_Y * nVoxel_Z)) * 100 / 2;
            msg = sprintf('%.2f percent. ', Prozent);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));

            for X = 1:nVoxel_X-1
                if (masksignal(Y, X, Z) ~= 0) && (masksignal(Y, X+1, Z) ~= 0)
                    for N1 = 1:nMinimaPerVoxel(Y, X, Z)
                        % Edges from a to t
                        [a1, a2, a3] = deal(0);
                        a1 = indexLocalMinimaRescale(Y, X, Z, N1);
                        if N1 ~= 1
                            a2 = indexLocalMinimaRescale(Y, X, Z, N1-1);
                        end
                        a3 = indexLocalMinimaRescale(Y, X+1, Z, nMinimaPerVoxel(Y, X+1, Z));
                        weight = g2(a1, a2, a3, @f);
                        node1 = NodeIndicies(Y, X, Z, N1);

                        if weight > 0 
                            InterColumnEdges(rem, 1) = node1;
                            InterColumnEdges(rem, 2) = T_Node_ID;
                            InterColumnEdges(rem, 3) = weight;
                            rem = rem + 1;
                            weight = 0;
                        end

                        for N2 = 2:nMinimaPerVoxel(Y, X+1, Z)
                            % Edges from a to b
                            [a1, a2, a3, a4] = deal(0);
                            a1 = indexLocalMinimaRescale(Y, X, Z, N1);
                            if N1 ~= 1
                                a2 = indexLocalMinimaRescale(Y, X, Z, N1-1);
                            end
                            a3 = indexLocalMinimaRescale(Y, X+1, Z, N2);
                            a4 = indexLocalMinimaRescale(Y, X+1, Z, N2-1);
                            weight = g(a1, a2, a3, a4, @f);
                            node1 = NodeIndicies(Y, X, Z, N1);
                            node2 = NodeIndicies(Y, X+1, Z, N2);


                            if (weight > 0) && (a3 ~= 0)
                                InterColumnEdges(rem, 1) = node1;
                                InterColumnEdges(rem, 2) = node2;
                                InterColumnEdges(rem, 3) = weight;
                                rem = rem + 1;
                                weight = 0;
                            end
                        end
                    end
                end
            end
        end
    end

    for Z = 1:nVoxel_Z
        for Y = 1:nVoxel_Y

            Prozent = 50 + ((((Z - 1) * nVoxel_Y) + Y) / (nVoxel_Y * nVoxel_Z)) * 100 / 2;
            msg = sprintf('%.2f percent. ', Prozent);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));

            for X = 1:nVoxel_X-1
                if (masksignal(Y, X, Z) ~= 0) && (masksignal(Y, X+1, Z) ~= 0)
                    for N1 = 1:nMinimaPerVoxel(Y, X+1, Z)

                        % Edges from b to t
                        [a1, a2, a3] = deal(0);
                        a1 = indexLocalMinimaRescale(Y, X+1, Z, N1);
                        if N1 ~= 1
                            a2 = indexLocalMinimaRescale(Y, X+1, Z, N1-1);
                        end
                        a3 = indexLocalMinimaRescale(Y, X, Z, nMinimaPerVoxel(Y, X, Z));
                        weight = g2(a1, a2, a3, @f);
                        node1 = NodeIndicies(Y, X+1, Z, N1);

                        if weight > 0
                            InterColumnEdges(rem, 1) = node1;
                            InterColumnEdges(rem, 2) = T_Node_ID;
                            InterColumnEdges(rem, 3) = weight;
                            rem = rem + 1;
                            weight = 0;
                        end  


                        for N2 = 2:nMinimaPerVoxel(Y, X, Z)
                            % Edges from b to a
                            [a1, a2, a3, a4] = deal(0);
                            a1 = indexLocalMinimaRescale(Y, X+1, Z, N1);
                            if N1 ~= 1
                                a2 = indexLocalMinimaRescale(Y, X+1, Z, N1-1);
                            end
                            a3 = indexLocalMinimaRescale(Y, X, Z, N2);
                            a4 = indexLocalMinimaRescale(Y, X, Z, N2-1);
                            weight = g(a1, a2, a3, a4, @f);

                            node1 = NodeIndicies(Y, X+1, Z, N1);
                            node2 = NodeIndicies(Y, X, Z, N2);

                            if (weight > 0) && (a3 ~= 0)
                                InterColumnEdges(rem, 1) = node1;
                                InterColumnEdges(rem, 2) = node2;
                                InterColumnEdges(rem, 3) = weight;
                                rem = rem + 1;
                                weight = 0;
                            end
                        end
                    end
                end
            end
        end
    end


    t1 = toc(t1);
    fprintf('Done! (%.2fs)\n', t1);


    %% Edges in Y-Direction
    fprintf('Calculating Edges in Y-Direction... ');
    t1 = tic;
    reverseStr = '';
    for Z = 1:nVoxel_Z
        for Y = 1:nVoxel_Y-1

            Prozent = ((((Z - 1) * (nVoxel_Y - 1)) + Y) / ((nVoxel_Y - 1) * nVoxel_Z)) * 100 / 2;
            msg = sprintf('%.2f percent. ', Prozent);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));

            for X = 1:nVoxel_X
                if (masksignal(Y, X, Z) ~= 0) && (masksignal(Y+1, X, Z) ~= 0)
                    for N1 = 1:nMinimaPerVoxel(Y, X, Z)
                        % Edges from a to t
                        [a1, a2, a3] = deal(0);
                        a1 = indexLocalMinimaRescale(Y, X , Z, N1);
                        if N1 ~= 1
                            a2 = indexLocalMinimaRescale(Y, X, Z, N1-1);
                        end
                        a3 = indexLocalMinimaRescale(Y+1, X, Z, nMinimaPerVoxel(Y+1, X, Z));
                        weight = g2(a1, a2, a3, @f);
                        node1 = NodeIndicies(Y, X, Z, N1);

                        if weight > 0
                            InterColumnEdges(rem, 1) = node1;
                            InterColumnEdges(rem, 2) = T_Node_ID;
                            InterColumnEdges(rem, 3) = weight;
                            rem = rem + 1;
                            weight = 0;
                        end 

                        for N2 = 2:nMinimaPerVoxel(Y+1, X, Z)
                            % forwards
                            [a1, a2, a3, a4] = deal(0);
                            a1 = indexLocalMinimaRescale(Y, X, Z, N1);
                            if N1 ~= 1
                                a2 = indexLocalMinimaRescale(Y, X, Z, N1-1);
                            end
                            a3 = indexLocalMinimaRescale(Y+1, X, Z, N2);
                            a4 = indexLocalMinimaRescale(Y+1, X, Z, N2-1);
                            weight = g(a1, a2, a3, a4, @f);
                            node1 = NodeIndicies(Y, X, Z, N1);
                            node2 = NodeIndicies(Y+1, X, Z, N2);

                            if (weight > 0) && (a3 ~= 0)
                                InterColumnEdges(rem, 1) = node1;
                                InterColumnEdges(rem, 2) = node2;
                                InterColumnEdges(rem, 3) = weight;
                                rem = rem + 1;
                                weight = 0;
                            end
                        end
                    end
                end
            end
        end
    end


    for Z = 1:nVoxel_Z
        for Y = 1:nVoxel_Y-1

            Prozent = 50 + ((((Z - 1) * (nVoxel_Y - 1)) + Y) / ((nVoxel_Y - 1) * nVoxel_Z)) * 100 / 2;
            msg = sprintf('%.2f percent. ', Prozent);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));

            for X = 1:nVoxel_X
                if (masksignal(Y, X, Z) ~= 0) && (masksignal(Y+1, X, Z) ~= 0)
                    for N1 = 1:nMinimaPerVoxel(Y+1, X, Z)

                        % Edges from b to t
                        [a1, a2, a3] = deal(0);
                        a1 = indexLocalMinimaRescale(Y+1, X, Z, N1);
                        if N1 ~= 1
                            a2 = indexLocalMinimaRescale(Y+1, X, Z, N1-1);
                        end
                        a3 = indexLocalMinimaRescale(Y, X, Z, nMinimaPerVoxel(Y, X, Z));
                        weight = g2(a1, a2, a3, @f);
                        node1 = NodeIndicies(Y+1, X, Z, N1);

                        if weight > 0
                            InterColumnEdges(rem, 1) = node1;
                            InterColumnEdges(rem, 2) = T_Node_ID;
                            InterColumnEdges(rem, 3) = weight;
                            rem = rem + 1;
                            weight = 0;
                        end 

                        for N2 = 2:nMinimaPerVoxel(Y, X, Z)
                            % backwards
                            [a1, a2, a3, a4] = deal(0);
                            a1 = indexLocalMinimaRescale(Y+1, X, Z, N1);
                            if N1 ~= 1
                                a2 = indexLocalMinimaRescale(Y+1, X, Z, N1-1);
                            end
                            a3 = indexLocalMinimaRescale(Y, X, Z, N2);
                            a4 = indexLocalMinimaRescale(Y, X, Z, N2-1);
                            weight = g(a1, a2, a3, a4, @f);

                            node1 = NodeIndicies(Y+1, X, Z, N1);
                            node2 = NodeIndicies(Y, X, Z, N2);

                            if (weight > 0) && (a3 ~= 0)
                                InterColumnEdges(rem, 1) = node1;
                                InterColumnEdges(rem, 2) = node2;
                                InterColumnEdges(rem, 3) = weight;
                                rem = rem + 1;
                                weight = 0;
                            end
                        end
                    end
                end
            end
        end
    end


    t1 = toc(t1);
    fprintf('Done! (%.2fs)\n', t1);

    %% Edges in Z-Direction

    if nVoxel_Z > 1

        fprintf('Calculating Edges in Z-Direction... ');
        t1 = tic;
        reverseStr = '';

        for Z = 1:nVoxel_Z-1
            for Y = 1:nVoxel_Y

            Prozent = ((((Z - 1) * nVoxel_Y) + Y) / (nVoxel_Y * (nVoxel_Z - 1))) * 100 / 2;
            msg = sprintf('%.2f percent. ', Prozent);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));

                for X = 1:nVoxel_X
                    if (masksignal(Y, X, Z) ~= 0) && (masksignal(Y, X, Z+1) ~= 0)
                        for N1 = 1:nMinimaPerVoxel(Y, X, Z)
                            % Edges from a to t
                            [a1, a2, a3] = deal(0);
                            a1 = indexLocalMinimaRescale(Y, X , Z, N1);
                            if N1 ~= 1
                                a2 = indexLocalMinimaRescale(Y, X, Z, N1-1);
                            end
                            a3 = indexLocalMinimaRescale(Y, X, Z+1, nMinimaPerVoxel(Y, X, Z+1));
                            weight = g2(a1, a2, a3, @f);
                            node1 = NodeIndicies(Y, X, Z, N1);

                            if weight > 0
                                InterColumnEdges(rem, 1) = node1;
                                InterColumnEdges(rem, 2) = T_Node_ID;
                                InterColumnEdges(rem, 3) = weight;
                                rem = rem + 1;
                                weight = 0;
                            end 

                            for N2 = 2:nMinimaPerVoxel(Y, X, Z+1)
                                % forwards
                                [a1, a2, a3, a4] = deal(0);
                                a1 = indexLocalMinimaRescale(Y, X, Z, N1);
                                if N1 ~= 1
                                    a2 = indexLocalMinimaRescale(Y, X, Z, N1-1);
                                end
                                a3 = indexLocalMinimaRescale(Y, X, Z+1, N2);
                                a4 = indexLocalMinimaRescale(Y, X, Z+1, N2-1);
                                weight = g(a1, a2, a3, a4, @f);
                                node1 = NodeIndicies(Y, X, Z, N1);
                                node2 = NodeIndicies(Y, X, Z+1, N2);

                                if (weight > 0) && (a3 ~= 0)
                                    InterColumnEdges(rem, 1) = node1;
                                    InterColumnEdges(rem, 2) = node2;
                                    InterColumnEdges(rem, 3) = weight;
                                    rem = rem + 1;
                                    weight = 0;
                                end
                            end
                        end
                    end
                end
            end
        end


        for Z = 1:nVoxel_Z-1
            for Y = 1:nVoxel_Y

            Prozent = 50 + ((((Z - 1) * nVoxel_Y) + Y) / (nVoxel_Y * (nVoxel_Z - 1))) * 100 / 2;
            msg = sprintf('%.2f percent. ', Prozent);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));

                for X = 1:nVoxel_X
                    if (masksignal(Y, X, Z) ~= 0) && (masksignal(Y, X, Z+1) ~= 0)
                        for N1 = 1:nMinimaPerVoxel(Y, X, Z+1)
                            % Edges from b to t
                            [a1, a2, a3] = deal(0);
                            a1 = indexLocalMinimaRescale(Y, X, Z+1, N1);
                            if N1 ~= 1
                                a2 = indexLocalMinimaRescale(Y, X, Z+1, N1-1);
                            end
                            a3 = indexLocalMinimaRescale(Y, X, Z, nMinimaPerVoxel(Y, X, Z));
                            weight = g2(a1, a2, a3, @f);
                            node1 = NodeIndicies(Y, X, Z+1, N1);

                            if weight > 0
                                InterColumnEdges(rem, 1) = node1;
                                InterColumnEdges(rem, 2) = T_Node_ID;
                                InterColumnEdges(rem, 3) = weight;
                                rem = rem + 1;
                                weight = 0;
                            end 

                            for N2 = 2:nMinimaPerVoxel(Y, X, Z)
                                % backwards
                                [a1, a2, a3, a4] = deal(0);
                                a1 = indexLocalMinimaRescale(Y, X, Z+1, N1);
                                if N1 ~= 1
                                    a2 = indexLocalMinimaRescale(Y, X, Z+1, N1-1);
                                end
                                a3 = indexLocalMinimaRescale(Y, X, Z, N2);
                                a4 = indexLocalMinimaRescale(Y, X, Z, N2-1);
                                weight = g(a1, a2, a3, a4, @f);

                                node1 = NodeIndicies(Y, X, Z+1, N1);
                                node2 = NodeIndicies(Y, X, Z, N2);

                                if (weight > 0) && (a3 ~= 0)
                                    InterColumnEdges(rem, 1) = node1;
                                    InterColumnEdges(rem, 2) = node2;
                                    InterColumnEdges(rem, 3) = weight;
                                    rem = rem + 1;
                                    weight = 0;
                                end
                            end
                        end
                    end
                end
            end
        end

    t1 = toc(t1);
    fprintf('Done! (%.2fs)\n', t1);

    end
    InterColumnEdges = InterColumnEdges(1:rem-1, :);

    %% Generate Adjacency Matrix from the Edgelists
    t = tic;
    fprintf('Creating Graph... ');
    A = sparse( InterColumnEdges(:, 1), InterColumnEdges(:, 2), InterColumnEdges(:, 3), nNodes+2, nNodes+2 );
    B = sparse( IntraColumnEdges(:, 1), IntraColumnEdges(:, 2), IntraColumnEdges(:, 3), nNodes+2, nNodes+2 );
    C = sparse( IntraColumnEdges(:, 2), IntraColumnEdges(:, 1), IntraColumnEdges(:, 4), nNodes+2, nNodes+2 );
    clearvars IntraColumnEdges InterColumnEdges
    D = A + B + C;
    clearvars A B C

    %% Create MATLAB-Graph and calculate maxflow
    G = digraph(D);
    clear D
    t = toc(t);
    fprintf('Done! (%.2fs)\n', t);
    fprintf('Solving Graph... ');
    
    t = tic;
    [~, ~, cs, ~] = maxflow(G, S_Node_ID, T_Node_ID);
    t = toc(t);
    
    fprintf('Done! (%.2fs)\n', t);

    %% Reconstruct fm
    fprintf('Reconstructing Fieldmap...');
    CS_Array = zeros(nVoxel_Y, nVoxel_X, nVoxel_Z, nMaxNodesPerVoxel);
    reverseStr = '';

    for Z = 1:nVoxel_Z
        for Y = 1:nVoxel_Y

            Prozent = ((((Z - 1) * nVoxel_Y) + Y) / (nVoxel_Y * nVoxel_Z)) * 100;
            msg = sprintf('%.2f percent. ', Prozent);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));

            for X = 1:nVoxel_X
                for N = 1:nMinimaPerVoxel(Y, X, Z)
                    if (masksignal(Y, X, Z) ~= 0)
                        if ismembc(NodeIndicies(Y, X, Z, N), cs)
                            CS_Array(Y, X, Z, N) = NodeIndicies(Y, X, Z, N);
                        end
                    end
                end
            end
        end
    end

    [~, IndexMax] = max(CS_Array, [], 4);

    fm = zeros(nVoxel_Y, nVoxel_X, nVoxel_Z);
    for Z = 1:nVoxel_Z
        for Y = 1:nVoxel_Y
            for X = 1:nVoxel_X
                if masksignal(Y, X, Z) ~= 0
                    fm(Y, X, Z) = indexLocalMinimaRescale(Y, X, Z, IndexMax(Y, X, Z));
                else
                    fm(Y, X, Z) = -algoParams.range_fm(1);
                end
            end
        end
    end
    fprintf('Done!\n');

    %% Backproject to physical Fieldmap
    dfm = algoParams.range_fm(1) + fm;

    %% Calculate r2starmap and water-fat separated images
    fprintf('Calculating r2starmap and water-fat separated images... ');
    water = ones(nVoxel_Y, nVoxel_X, nVoxel_Z);
    fat = ones(nVoxel_Y, nVoxel_X, nVoxel_Z);
    r2starmap = ones(nVoxel_Y, nVoxel_X, nVoxel_Z);  
    for Z = 1:nVoxel_Z
        
        msg = sprintf('%.2f percent. ', Z / nVoxel_Z * 100);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));

        tmp_imDataParams = imDataParams;
        tmp_imDataParams.images = tmp_imDataParams.images(:, :, Z, :, :);
        
        r2starmap(:, :, Z) = estimateR2starGivenFieldmap( tmp_imDataParams, algoParams, squeeze(dfm(:, :, Z)) );
        amps = decomposeGivenFieldMapAndDampings( tmp_imDataParams, algoParams, squeeze(dfm(:, :, Z)), r2starmap(:, :, Z), r2starmap(:, :, Z) );
        
        waterimage = squeeze(amps(:, :, 1, :));
        fatimage = squeeze(amps(:, :, 2, :));

        water(:, :, Z) = waterimage;
        fat(:, :, Z) = fatimage;
    end
    fprintf('Done!\n');

    %% Transfer results into outParams structure
    outParams.water = water .* masksignal;
    outParams.fat = fat .* masksignal;
    outParams.fieldmap = dfm;
    outParams.r2starmap = r2starmap .* masksignal;
end