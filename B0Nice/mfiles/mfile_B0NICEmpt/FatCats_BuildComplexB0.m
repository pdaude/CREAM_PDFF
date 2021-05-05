%
function [algoParams] = FatCats_BuildComplexB0(algoParams)
%
complex_image = algoParams.complex_image;
complex_image(isnan(complex_image)) = 1;
matrix_size = algoParams.matrix_size;
flag_B0 = algoParams.index_B0(1);
index_start = algoParams.index_B0(2);
index_end = algoParams.index_B0(3);
%
complex_B0 = zeros(matrix_size(1),matrix_size(2),matrix_size(3));
    %
    if flag_B0 == 1  
        complex_B0(:,:,:) = complex_image(:,:,:,1,index_end).*conj(complex_image(:,:,:,1,index_start));
    elseif flag_B0 == 2 
        complex_0 = complex_image(:,:,:,1,index_start);
        complex_1 = complex_image(:,:,:,1,index_start+1); 
        complex_2 = complex_image(:,:,:,1,index_end);
        complex_10 = complex_1.*conj(complex_2);
        complex_20 = complex_0.*conj(complex_2);
        complex_B0 = conj(complex_10.*complex_20);
    %elseif flag_B0 ==3
    %    complex_0 = complex_image(:,:,:,1,index_start);
    %    complex_1 = complex_image(:,:,:,1,index_end);   
    %    complex_tmp = complex_1.*conj(complex_0);
    %    complex_tmp = exp(1i*angle(complex_tmp));
    %    complex_tmp = complex_tmp.^2;
    %    phase_tmp = angle(complex_tmp);
    %    complex_B0 = abs(complex_1).*exp(1i*phase_tmp); 
    end
    %
    if flag_B0 ~=3
        phase_B0 = angle(complex_B0);
        if abs(algoParams.phase_echo(1)) <= abs(algoParams.phase_echo(2))
        complex_B0 = exp(1i*phase_B0).*abs(complex_image(:,:,:,1,1));
        else
        complex_B0 = exp(1i*phase_B0).*abs(complex_image(:,:,:,1,2));    
        end
    end
    %
    if algoParams.FILTsize > 1
        H = fspecial('average',[algoParams.FILTsize algoParams.FILTsize]);
        for index_slice = 1:algoParams.matrix_size(3)
        complex_B0slice(:,:) = complex_B0(:,:,index_slice);
        real_B0slice = imfilter(real(complex_B0slice),H);
        imag_B0slice = imfilter(imag(complex_B0slice),H);
        complex_B0(:,:,index_slice) = real_B0slice + 1i*imag_B0slice;
        end
    end
    algoParams.complex_B0 = complex_B0;
        
