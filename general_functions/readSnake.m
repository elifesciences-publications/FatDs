function snake = readSnake(snakeFile, nTimePts)

    fid = fopen(snakeFile);
    
    if fid == -1
        error(['file not found: ' snakeFile]);
    end

    % scan through the header
    tline = fgets(fid);
    while ischar(tline) && ~strcmp(tline(1), '#')
        tline = fgets(fid);
    end
    tline = fgets(fid); % skip the 0

    % read the snake
    % lines are of the form [t p x y 0]
    % t is frame number, p contour index
    snake = cell([nTimePts 1]);
    n = 1;
    
    tline = fgets(fid);
    while ischar(tline) 

        if ~strcmp(tline(1), '#')
            numline = str2num(tline);
            t = uint16(numline(1));
            snake{t, n} = cat(1,snake{t, n}, numline(:,3:4) + 1);
            tline = fgets(fid);
            
        % start a new snake
        else
            n = n + 1;
            snake = [snake cell([nTimePts 1])];
            fgets(fid); % skip the #
            tline = fgets(fid); % skip the 0
        end
    end
    
end