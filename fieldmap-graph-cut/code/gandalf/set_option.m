function val = set_option(Options, fieldName, defaultVal)
    
    if nargin < 3
        defaultVal = [];
    end
    
    if isfield(Options, fieldName) & ~isempty(Options.(fieldName))
        val = Options.(fieldName);

    else
        val = defaultVal;
    end