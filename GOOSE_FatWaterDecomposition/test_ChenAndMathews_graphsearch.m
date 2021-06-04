%% test_ChenAndMathews
%%
%% Test fat-water algorithms on (ISMRM Challege's) data
%%
%% Author: Chen Cui, Mathews Jacob
%% Date created: March 11, 2013
%% Date last modified: March 15, 2013

clear all
%% Add to matlab path
% BASEPATH = './';
% addpath([BASEPATH 'common/']);

%% start the algorithm
dataname = '01.mat';
dataname ='/home/pdaude/Projet_Python/CREAM_PDFF/fw_i3cm1i_3pluspoint_berglund_QPBO/test_cases/01.mat'
z = 1; % processing the zth slice -- the third dimension of original images data
tic
[outParams] = GraphSearch(dataname,z);
toc