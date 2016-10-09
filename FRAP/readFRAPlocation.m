function [x y w h] = readFRAPlocation(filepath, rgnfname)
    % read FRAP region from .rgn file
    %
    % [x y w h] = readFRAPlocation(filepath, rgnfname)
    % 
    % x,y:  location of upper left corner
    % w,h:  width, height of box?

    % read last line of file, in clunky way
    fid = fopen(fullfile(filepath, rgnfname));
    tline = fgets(fid);
    prevline = [];
    while ischar(tline)
        tline = fgets(fid);
        if tline ~= -1 
            prevline = tline;
        end
    end
    
    % find position of equals sign in string
    N = regexp(prevline,'\=');
    % read out the part after the equals sign and convert to matrix
    X = str2num(prevline(23:end));

    x = X(1,1);
    y = X(1,2);
    w = X(2,1)-X(1,1);
    h = X(2,2)-X(1,2);

end