function echoMIP = get_echoMIP(signal)
% echoMIP = get_echoMIP(signal)

    magnitude = abs(signal);
    echoMIP = squeeze(sqrt(sum(magnitude.^2, ndims(signal))));

end
        