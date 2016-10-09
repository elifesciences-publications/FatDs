function [c, cx, cy, nx, ny] = broadImprofile(im, xi, yi, t)
    % like improfile but with lines that are more than a pixel thick and
    % without giving a number of pts
    %
    % [c, cx, cy] = fatImprofile(im, xi, yi, t)
    %
    % t:        line thickness, odd integer
    %
    % Idse Heemskerk, August 2014

    if nargin < 4
        t = 1;
    end
    
    if rem(t,2) == 0
        t = t+1;
        disp(['making t odd: ' num2str(t)]);
    end
    
    if size(xi, 2)==1
        xi = xi';
        yi = yi';
    end

    dy = gradient(yi);
    dx = gradient(xi);
    
    N = sqrt(dx.^2 + dy.^2);
    nx = -dy./N; % used to say dy/N (before Sep 2016)
    ny = dx./N;
    
    c = zeros([t length(xi)]);
    cx = zeros([t length(xi)]);
    cy = zeros([t length(xi)]);

    for i = -(t-1)/2:(t-1)/2
        
        idx = i + (t-1)/2 + 1;
        cx(idx,:) = xi + i*nx;
        cy(idx,:) = yi + i*ny;
        c(idx,:) = interp2(double(im), cx(idx,:), cy(idx,:), 'linear');
    end
end