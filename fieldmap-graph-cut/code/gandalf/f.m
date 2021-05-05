function value = f(r1,r2)
    value = 0;
    if r1 > r2
        value = abs(r1-r2)*abs(r1-r2);
    end
end