function VARPROparamsout = checkParamsAndSetDefaults_GANDALF(imDataParams, algoParams)

    validParams = 1;

    % Start by checking validity of provided data and recon parameters
    if size(imDataParams,3) > 1
      disp('ERROR: 2D recon -- please format input data as array of size SX x SY x 1 X nCoils X nTE')
      validParams = 0;
    end

    if length(algoParams.species) > 2
      disp('ERROR: Water=fat recon -- use a multi-species function to separate more than 2 chemical species')
      validParams = 0;
    end

    if length(imDataParams.TE) < 3
      disp('ERROR: 3+ point recon -- please use a different recon for acquisitions with fewer than 3 TEs')
      validParams = 0;
    end
    VARPROparamsout.validParams = validParams;
    
    dTE = diff(imDataParams.TE);
    if sum(abs(dTE - dTE(1)))<1e-6 % If we have uniform TE spacing
        dt = imDataParams.TE(2) - imDataParams.TE(1);
        fprintf('Uniform TE spacing: Period = TE(2) - TE(1)\n');
    else
        dt = imDataParams.TE(4) - imDataParams.TE(3);
        fprintf('Non-uniform TE spacing detected, assuming 2 UTE echos: Period = TE(4) - TE(3)\n');
    end    
    
    VARPROparamsout.species = algoParams.species;
    VARPROparamsout.useCUDA = set_option(algoParams, 'useCUDA', 1);
    VARPROparamsout.range_r2star = set_option(algoParams, 'range_r2star', [0 500]);
    VARPROparamsout.NUM_R2STARS = set_option(algoParams, 'NUM_R2STARS', 26);
    VARPROparamsout.sampling_stepsize = set_option(algoParams, 'sampling_stepsize', 2);
    VARPROparamsout.nSamplingPeriods = set_option(algoParams, 'nSamplingPeriods', 1);
    VARPROparamsout.airSignalThreshold_percent = set_option(algoParams, 'airSignalThreshold_percent', 5);
    VARPROparamsout.period = abs(1/dt);
    
    if isfield(algoParams, 'range_fm')
        VARPROparamsout.range_fm = algoParams.range_fm;
        VARPROparamsout.nSamplingPeriods = abs(diff(VARPROparamsout.range_fm)) / VARPROparamsout.period;
    else
        VARPROparamsout.nSamplingPeriods = set_option(algoParams, 'nSamplingPeriods', 1);
        VARPROparamsout.range_fm = [(-VARPROparamsout.nSamplingPeriods * VARPROparamsout.period / 2) (2*VARPROparamsout.sampling_stepsize + VARPROparamsout.nSamplingPeriods * VARPROparamsout.period / 2)];
    end
    
    disctreizationintervall = ceil( diff(VARPROparamsout.range_fm) );
    Numlayers = ceil( disctreizationintervall / VARPROparamsout.sampling_stepsize );
    VARPROparamsout.NUM_FMS = Numlayers;
    t = linspace(VARPROparamsout.range_fm(1), VARPROparamsout.range_fm(2), VARPROparamsout.NUM_FMS);
    gridspacing = t(2)-t(1);
    VARPROparamsout.gridspacing = gridspacing;
end