function [fatFrac_3D] = FatCats_getSpecies(algoParams)

%% Description 
% This function calculates the water, fat components and fat fraction map
% given the field map and R2_star
%
% % Inputs:
%           field:  Field map [rows, cols]
%           R2s:    R2_star map [rows, cols]
%           sig:    Acquired signal [rows, cols, number of echoes]
%           TE:     Echo times (sec) [number of echoes, 1]
% % Outputs:
%           water, fat and fat fraction components
%
% % % % % % % % % % % % % % % % % % % % % 
%       Abraam Soliman (asoliman@uwo.ca)
%       Western University
%       Feb 15, 2013
% % % % % % % % % % % % % % % % % % % % %
% modified by Junmin liu Feb 15, 2013

%% Calculate species (RHO)
pi_2 = 2*pi;
field_3D = algoParams.B0map_Hz;
R2star_map = algoParams.R2star_map_raw;
complex_image = algoParams.complex_image;
TE = algoParams.TE_seq;
error_bias = algoParams.error_bias;
amp_all = algoParams.model_r;
freq_all = algoParams.model_f;
freq_all = freq_all - error_bias;
%
if length(TE) > algoParams.CutOff
complex_image(:,:,:,:,(algoParams.CutOff+1):length(TE)) = [];
TE((algoParams.CutOff+1):length(TE)) = [];
end
%
matrix_size = size(complex_image);
for index_slice = 1:matrix_size(3)
field(:,:) = field_3D(:,:,index_slice); 
field = medfilt2(field,[5 5]);
R2s(:,:) = R2star_map(:,:,index_slice);
R2s = medfilt2(R2s,[5 5]);
sig(:,:,:) = complex_image(:,:,index_slice,1,:); 
%            

[rows cols] = size(field);
Rho = zeros(rows, cols, 2);
complexPhi = field + 1i*R2s./(2*pi);
matA = getMatA(TE, freq_all,amp_all);
Ainv = (matA'*matA) \ matA';

for j = 1:rows
    for i = 1:cols
        Psi_inv = diag(exp(-1i*2*pi* complexPhi(j, i).*TE));
        Rho(j, i, :) = Ainv * Psi_inv * squeeze(sig(j, i, :));
    end
end


%% Calculate fat fraction
% % Without noise-bias correction

water = Rho(:, :, 1);
fat = Rho(:, :, 2);
fatFrac = abs(fat)./(abs(fat) + abs(water));

%showWat = imadjust(abs(water) ./ max(max(abs(water))));
%showFat = imadjust(abs(fat) ./ max(max(abs(fat))));
 
%figure, imagesc(showWat), title('Water'), colormap(gray)
%figure, imagesc(showFat), title('Fat'), colormap(gray)
%figure, imagesc(fatFrac, [0 1]), title('Fat Fraction'), colormap(gray)
%
fatFrac_3D(:,:,index_slice) = fatFrac;
%index_slice
%pause
end
%% Function: get A_Matrix (6 peaks model)

function matA = getMatA(t, freq_all,amp_all)

matA = ones(length(t), 2);

for k = 1:length(t)
    matA(k, 2) = exp(1i* 2*pi* freq_all.* t(k)) * amp_all';
end

