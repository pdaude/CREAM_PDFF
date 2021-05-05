function current = ICM(prev, nB0, maxICMupdate, nICMiter, J, V, wx, wy, wz, left, right, up, down, above, below, lambda)
%ICM Iterated conditional modes
%   Performs iterated conditional modes

current = prev;
for k = 1:nICMiter  % ICM iterate
    prev = current;
    min_cost = Inf(size(current));
    
    updates = zeros(1,2*maxICMupdate+1);  % Update order
    % Odd are positive
    updates(3:2:end) = 1:maxICMupdate;
    % Even are negative
    updates(2:2:end) = -(1:maxICMupdate);
    for update = updates
        % Unary cost:
        cost = lambda * J(sub2ind(size(J), mod((prev+update)-1, nB0)+1, (1:size(J,2)).'));
        % Binary costs:
        cost(right) = cost(right) + wx.*V(abs(mod((prev(right)+update-1), nB0)+1-prev(left))+1);
        cost(left) = cost(left) + wx.*V(abs(mod((prev(left)+update-1), nB0)+1-prev(right))+1);
        cost(down) = cost(down) + wy.*V(abs(mod((prev(down)+update-1), nB0)+1-prev(up))+1);
        cost(up) = cost(up) + wy.*V(abs(mod((prev(up)+update-1), nB0)+1-prev(down))+1);
        if ~isempty(wz)
            cost(below) = cost(below) + wz.*V(abs(mod((prev(below)+update-1), nB0)+1-prev(above))+1);
            cost(above) = cost(above) + wz.*V(abs(mod((prev(above)+update-1), nB0)+1-prev(below))+1);
        end
        
        current(cost < min_cost) = mod(prev(cost < min_cost)+update-1, nB0)+1;
        min_cost(cost < min_cost) = cost(cost < min_cost);
    end
    
end

end