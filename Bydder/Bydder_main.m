function [params, sse,sse1,ym] = Bydder_main(imDataParams,algoParams)
X = fieldnames(algoParams);
Y = struct2cell(algoParams);
C = [X,Y].';
[nx,ny,nz,ncoils,ne]=size(imDataParams.images);
data = reshape(imDataParams.images,nx,ny,nz,ne);
[params,sse1]=presco(imDataParams.TE,data,imDataParams.FieldStrength, C{:});
[sse,ym]=calculate_residual(params,algoParams,imDataParams);
end