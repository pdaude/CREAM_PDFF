function minima = findMinima(f)
%FINDMINIMA Returns all local minima
%   Returns all local minima of the 'periodic' vector f

minima = find(f<[f(end);f(1:end-1)] & f<[f(2:end);f(1)]);

end