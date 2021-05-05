function tissueMask = get_tissueMask(signal, airSignalThreshold_percent)
% tissueMask = get_tissueMask(signal, airSignalThreshold_percent)

    echoMIP = get_echoMIP(signal);

    threshold = airSignalThreshold_percent/100 * max(echoMIP(:));

    tissueMask = echoMIP >= threshold .* ones(size(echoMIP));

end