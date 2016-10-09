function newSnake = growSnake(snake, L)
    % grow snake by L pixels on both ends

    dx = gradient(snake(:,1));
    dy = gradient(snake(:,2));
    N = sqrt(dx.^2 + dy.^2);

    % proper distance / positions to extend
    s = 1:L;

    % number of positions to use in approximating tangent
    n = 2;
    
    dxavg = mean(dx(end-n:end)./N(end-n:end));
    dyavg = mean(dy(end-n:end)./N(end-n:end));

    newEnd = [snake(end,1) + s'*dxavg, snake(end,2) + s'*dyavg];

    dxavg = mean(dx(1:1+n)./N(1:1+n));
    dyavg = mean(dy(1:1+n)./N(1:1+n));

    newStart = [snake(1,1) - s'*dxavg, snake(1,2) - s'*dyavg];

    newSnake = cat(1, flipud(newStart), snake, newEnd);
end