% ReadMe

% Author: Chen Cui; Xiaodong Wu; John Newell; Mathews Jacob
% Date: 02/03/2014

% This toolbox is an implementation of paper "Fat water decomposition using 
% GlObally Optimal Surface Estimation (GOOSE) algorithm".

% It is designed as form of blackbox with:
% input: Matlab format, such as "01.mat" in 2012 ISMRM Challenge;
% output: "OutParams", including fielmap, separated water/fat images, and 
%         r2starmap in our case. (form in ISMRM W/F Workshop)

% How to run:
% 1. Open file "test_ChenAndMathews_graphsearch.m";
% 2. Change the dataname to the one that needs to be processed;
% 3. Specify which slice in z-direction will be passed into algorithm box;
%    z-direction is the 3rd dimension in "imDataParas.images".
% Then you are good to go!

% The structure of the algorithm is in function "GraphSearch.m", which
% contains five steps:
% 1. load the data;
% 2. Set reconstruction parameters;
% 3. Compute VARPRO residue;
% 4. Estimate fieldmap by graph search;
% 5. Retrieve fat and water images;

% The core part of the algorithm is done by
% "detect_optimal_surface.mexa64". It has been successfully run on 
% OS: Linux 3.4.11-2.16-desktop x86_64 (GLNXA64),
% System: openSUSE 12.2 (x86_64)
% CPU: Intel(R) Xeon(R) CPU W3565 @ 3.20GHz
% RAM:  23.6 GiB. 
% However, a miniumum of 8 GiB is enough for the toolbox to execute. 
% We also tested on Windows operating system after the Challenge
% competition, therefore, it should work on both Linux and Windows
% machines.

% We also attached the original "detect_optimal_surface.cpp" file and the
% graph library "optnet", just in case it has to be recompiled. 
  
