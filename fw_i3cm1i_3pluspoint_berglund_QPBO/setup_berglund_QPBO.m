%% setup fw_i3cm1i_3pluspoint_berglund_QPBO
%
% add parent folder of fw_i3cm1i_3pluspoint_berglund_QPBO.m to MATLAB path
%

%% parse location of setup_berglund_QPBO.m
[pathstr,name,ext] = fileparts(mfilename('fullpath'));

%% Store present working directory and change directory temporarily to matlabroot
pwd_path = pwd;
cd(matlabroot);

if exist('fw_i3cm1i_3pluspoint_berglund_QPBO')~=2,
    
    %% add path where setup and test scripts reside
    addpath(pathstr);
    
    %% add deeper path to fw_i3cm1i_3pluspoint_berglund_QPBO.m
    pathstr_to_add = sprintf('%s/berglund/QPBO',pathstr);
    addpath(pathstr_to_add);
    
    %% save path
    savepath;
    
    %% check again
    if exist('fw_i3cm1i_3pluspoint_berglund_QPBO')~=2,
        disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO.m appears to NOT be setup. Check that ''%s'' is added to MATLAB path.', pathstr_to_add) );
    else
        disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO.m is now accessible at ''%s''', which('fw_i3cm1i_3pluspoint_berglund_QPBO') ) );
    end
else
    disp( sprintf('fw_i3cm1i_3pluspoint_berglund_QPBO.m is already accessible at ''%s''', which('fw_i3cm1i_3pluspoint_berglund_QPBO') ) );
end

%% change back to the original directory
cd(pwd_path);
