function [ algoParams] = initial_b0MAPPING(algoParams )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
   [algoParams] = mask_generation4PUROR(algoParams); 
%calculating initial B0
    [unwrap_phase] = PUROR2D_B0Matlab(algoParams.complex_B0,algoParams.mask4unwrap,algoParams.mask4supp,algoParams.mask4STAS);       
%-----------------------------------------
%
    algoParams.B0map_Hz = unwrap_phase./(2*pi*algoParams.delta_TE4B0);

