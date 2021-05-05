function [arrayScaled, rescaleSlope, rescaleIntercept] = scale_array2interval(array, intervalVector)
% [arrayScaled, rescaleSlope, rescaleIntercept] = scale_array2interval(array, intervalVector)
% 
% AUTHOR:      Maximilian N. Diefenbach <maximilian.diefenbach@tum.de>
% AFFILIATION: Body Magnetic Resonance Research Group
%              Department of Diagnostic and Interventional Radiology
%              Technical University of Munich
%              Klinikum rechts der Isar
%              Ismaninger Str. 22, 81675 Muenchen
% URL:         http://www.bmrrgroup.de

    array = double(array);
    
    % scale by slope
    newValueRange = max(intervalVector(:)) - min(intervalVector(:));
    valueRange = max(array(:)) - min(array(:));
    if valueRange == 0
        rescaleSlope = 1;
    else
        rescaleSlope = newValueRange / valueRange;
    end
    arrayScaled = rescaleSlope * array;
    
    % shift by intercept
    rescaleIntercept = min(intervalVector(:)) - min(arrayScaled(:));
    arrayScaled = arrayScaled + rescaleIntercept;
    
    % test
    tol = 1e-3;
    assert(abs(min(arrayScaled(:)) - min(intervalVector(:))) <= tol);
    assert(abs(max(arrayScaled(:)) - max(intervalVector(:))) <= tol);

end