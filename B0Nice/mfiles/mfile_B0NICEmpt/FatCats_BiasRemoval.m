function [outParams] = FatCats_BiasRemoval( algoParams )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%#########################################  
    [FF4bias]= FatCats_getSpecies(algoParams);
    %
    outParams.fm_noFF = algoParams.B0map_Hz; 
    algoParams.B0map_Hz = algoParams.B0map_Hz - algoParams.error_bias.*FF4bias;
      
outParams.fm = algoParams.B0map_Hz;
outParams.FF = FF4bias;
outParams.R2star = algoParams.R2star_map_raw;
%