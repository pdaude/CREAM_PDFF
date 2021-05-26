%
function [algoParams] = FatCats_VerifyInputData(algoParams)
%
B0_strength = algoParams.B0_strength; % unit:T
spin_dir = algoParams.PrecessionIsClockwise;
TE_seq = algoParams.TE_seq;% unit:s
%
 
model_f = algoParams.species(2).frequency; % Should be in Hz 
model_r = algoParams.species(2).relAmps; % Should be Normalize
% model_r = [0.0870 0.6930 0.1280 0.0040 0.0390 0.0480];
% model_f = [-242.7060 -217.1580 -166.0620 -123.9078 -24.9093 38.3220];
% model_f = model_f.*(B0_strength/1.5);
%
algoParams.model_r = model_r;
algoParams.model_f = model_f;
%
%------------------------------
phase_echo = zeros(1,length(TE_seq));
mag_echo = zeros(1,length(TE_seq));
    for index_echo = 1:length(TE_seq)
    TE_tmp = TE_seq(index_echo);
    model_complex = model_r.*exp(1i.*(2*pi*TE_tmp.*model_f));
    model_complex_sum = sum(model_complex);
    phase_echo(index_echo) = angle(model_complex_sum);
    mag_echo(index_echo) = abs(model_complex_sum);
    end
%phase_echo/pi
%close all
%plot(phase_echo/pi,'o-r');
%pause
algoParams.Phase4Echo_ori = phase_echo;
algoParams.Mag4Echo_ori = mag_echo;
%-----------------------------
%
complex_3D = algoParams.complex_image;
matrix_size = size(complex_3D);
%
complex_image = zeros(matrix_size(1),matrix_size(2),matrix_size(3),1,matrix_size(5));
%
if matrix_size(4) > 1
    for index_slice = 1:matrix_size(3)
    complex_slice(:,:,1,:,:) = complex_3D(:,:,index_slice,:,:);
    [sliceComb] = coilCombine(complex_slice);
    im_tmp(:,:) = sliceComb(:,:,1,1,1);
        if index_slice == matrix_size(3)
        figure(1)
        subplot(1,2,1);imagesc(abs(im_tmp));
        subplot(1,2,2);imagesc(angle(im_tmp));
        index_slice
        pause
        end
    complex_image(:,:,index_slice,1,:) = sliceComb;
    end
else
    complex_image = complex_3D;
end
clear complex_3D
%
if spin_dir == -1
complex_image = conj(complex_image);
end
%
algoParams.complex_image = complex_image;
matrix_size = size(complex_image);
algoParams.matrix_size = matrix_size;

%