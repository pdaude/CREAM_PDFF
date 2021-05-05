%
function [algoParams] = FatCats_FFfromPhase(algoParams)
%
pi_2 = 2*pi;
%
complex_image = algoParams.complex_image;
B0_map = (pi_2*algoParams.delta_TE4B0).*algoParams.B0map_Hz;

matrix_size = size(complex_image);
%---------------------------------------------------------
    if algoParams.index_B0(1) == 3
        index_start = 1;
        echo_tmp = 2;
    elseif algoParams.index_B0(1) == 1 && (algoParams.index_B0(3) - algoParams.index_B0(2)) == 1
        index_start = algoParams.index_B0(3);
        echo_tmp = (algoParams.index_B0(3)+1);
    elseif algoParams.index_B0(1) == 2
        index_start = algoParams.index_B0(2);
        echo_tmp = (algoParams.index_B0(3));
    else
        index_start = algoParams.index_B0(2);
        echo_tmp = (algoParams.index_B0(3)-1);
    end
    %index_start
    %echo_tmp
    
    complex_B0Map_table = zeros(matrix_size(1),matrix_size(2),matrix_size(3),(matrix_size(5)-index_start));
for index_slice = 1:matrix_size(3)
    B0_unit(:,:) = B0_map(:,:,index_slice);
    jj = 1;
    %
    for index_echo = (index_start+1):echo_tmp
        TEratio_echo = (algoParams.TE_seq(index_echo)-algoParams.TE_seq(index_start))/algoParams.delta_TE4B0;
        if mod(TEratio_echo,0.5) > 0.1 || algoParams.index_B0(1) == 3 || (algoParams.index_B0(3) - algoParams.index_B0(2)) == 2
        Theoritical_FatPhase = angle(exp(1i*(algoParams.phase_echo(index_echo)))*exp(-1i*algoParams.phase_echo(index_start)));
        %Theoritical_FatPhase./pi;
        B0Map_tmp = B0_unit.*TEratio_echo;
        B0Map_tmp = medfilt2(B0Map_tmp,[3 3]);
        B0Map_tmp = exp(-1i*B0Map_tmp);
        complex_B0Map_table(:,:,index_slice,(index_echo - index_start)) = B0Map_tmp; 
        %
        complex_HP = complex_image(:,:,index_slice,index_echo).*conj(complex_image(:,:,index_slice,index_start));
        complex_HP = complex_HP.*B0Map_tmp;
        phaseFF_tmp(:,:,index_slice,jj) = angle(complex_HP)./Theoritical_FatPhase;
        %
        %close all
        %figure(2)
        %subplot(1,3,1);imagesc(angle(B0Map_tmp),[-pi pi]);colormap jet;axis square;axis off;
        %subplot(1,3,2);imagesc(angle(complex_HP),[-pi pi]);colormap jet;axis square;axis off;
        %subplot(1,3,3);imagesc(angle(complex_HP)./Theoritical_FatPhase,[-3 3]);colormap jet;axis square;axis off;
        %index_echo
        %index_slice
        %pause
        jj = jj+1;
        end
    end
end
%

    algoParams.phaseFF = mean(phaseFF_tmp,4);
    if algoParams.index_B0(1) == 3
    algoParams.phaseFF = abs(algoParams.phaseFF);
    end
end
%------------------------------------------------------
%figure(1);imshow3D(algoParams.phaseFF,[0 1]);colormap gray;
%pause
%



