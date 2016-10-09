function [data nChannels] = readCombineStack(filepath, fname) 
    % read the data from separate channel stacks
    % fname:    cell array with filenames of the separate channels

    % load the Bio-Formats library into the MATLAB environment
    autoloadBioFormats = 1;
    status = bfCheckJavaPath(autoloadBioFormats);
    assert(status, ['Missing Bio-Formats library. Either add loci_tools.jar '...
        'to the static Java path or add it to the Matlab path.']);

    % initialize logging
    loci.common.DebugTools.enableLogging('INFO');
    
    % knowing that we have multiple files each storing one channel:
    nChannels = numel(fname); %reader.getSizeC(); 

    % create bioformats reader for file
    reader = {};
    for c = 1:nChannels
        fullfname = fullfile(filepath, fname{c});
        disp(fullfname);
        reader{c} = bfGetReader(fullfname);
    end

    % dimensions (should match between files so just take first)
    nTimePts = reader{1}.getSizeT();
    xSize = reader{1}.getSizeX();
    ySize = reader{1}.getSizeY();
    zSize = reader{1}.getSizeZ();

    data = zeros([ySize xSize zSize 3 nTimePts], 'uint16');

    for t = 1:nTimePts
        for z = 1:zSize
            for c = 1:nChannels
                % convert multidimensional index to linear index
                % order is usually x y c z t, check
                % (perhaps include reader.isOrderCertain())
                if strcmp(reader{1}.getDimensionOrder(), 'XYCZT')
                    i = z + (t-1)*zSize;
                else
                    error('dimension order is not XYCZT, change script');
                end

                % read slice
                data(:,:,z,c,t) = bfGetPlane(reader{c}, i);
            end

            % progress indicator
            fprintf('.');
            if rem(i,80) == 0
                fprintf('\n');
            end
        end
    end
    fprintf('\n');

end
