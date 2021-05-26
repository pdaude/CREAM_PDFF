clear all 
clc
close all
%Add Paths :
%/home/pdaude/Projet_Python/CREAM_PDFF/Bydder with subfolders

algoParams=ReadYaml('/home/pdaude/Projet_Python/CREAM_PDFF/algoParams/Bydder_algoParams.yml');
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 4.7;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [0.9,1.3,2.1,2.76,4.31,5.3];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
algoParams.gyro=42.57747892;

X = fieldnames(algoParams);
Y = struct2cell(algoParams);
C = [X,Y].';


img=load('/home/pdaude/Projet_Python/CREAM_PDFF/fw_i3cm1i_3pluspoint_berglund_QPBO/test_cases/01.mat');
img.imDataParams.images = single(img.imDataParams.images);% Could not be integer
%doesn't take into consideration Precession is clockwise
[nx,ny,nz,ncoils,ne]=size(img.imDataParams.images);
data = reshape(img.imDataParams.images,nx,ny,nz,ne);
[params,sse]=presco(img.imDataParams.TE,data,img.imDataParams.FieldStrength, C{:});
[sse2,ym]=calculate_residual(params,algoParams,img.imDataParams);
save('/home/pdaude/Projet_Python/CREAM_PDFF/Bydder/output_Bydder_01.mat','params','sse','ym','sse2');
