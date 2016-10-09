function [fnames, lims, label] = preIlastik(datapath, dataset, corder, saveims, imposesLim)
    % save Fat-Ds and nuclei with same LUT for many folders
    %
    % [fnames, lims, label] = preIlastik(filepath, saveims)
    %
    % filepath: cell array of file paths
    % saveims:  boolean, actually write tifs to disk?
    %
    % fnames:   filename of raw (SIP) data in each folder
    % lims:     LUT limits 
    % label:    label for each image (based on SIP filename)

    fnames = {};
    lims = {};
    label = {};

    ims = {};
    
    segResultsDir = 'matlab_seg';
    
    % regenerate filepath and batchIdx
    filepath = {};
    nDatasets = 1;
    for i = 1:numel(dataset)
        for j = 1:numel(dataset{i})
            filepath{nDatasets} = fullfile(datapath{i}, dataset{i}{j});
            nDatasets = nDatasets + 1;
        end
    end
    
    n = 1;
    batchIdx = {};
    for i = 1:numel(datapath)
        batchIdx{i} = n:(n+numel(dataset{i})-1);
        n = n+numel(dataset{i});
    end
    
    % make directories for output at batch level
    for i = 1:numel(batchIdx)
        if ~exist(fullfile(datapath{i}, 'nuclei'), 'dir')
            mkdir(fullfile(datapath{i}, 'nuclei'));
        end
        if ~exist(fullfile(datapath{i}, 'FatDs'), 'dir')
            mkdir(fullfile(datapath{i}, 'FatDs'));
        end
    end
    
    % scan directories voor data
    for j = 1:numel(batchIdx)
        
        for i = batchIdx{j}
            
            found = false;

            % identify file starting with SUM
            %----------------------------

            dircont = dir(filepath{i});
            for j = 1:numel(dircont)

                dirj = dircont(j).name;

                if regexp(dirj,'AVG') == 1

                    fname = fullfile(filepath{i}, dircont(j).name);
                    fnames{i} = fname;
                    % THIS ASSUMES SOME CONVENTION 
                    label{i} = dirj(5:(regexp(dirj,'MM')-2));
                    found  = true;
                end
            end
            if ~found 
                for j = 1:numel(dircont)
                    dirj = dircont(j).name;

                    if regexp(dirj,'SUM') == 1

                        fname = fullfile(filepath{i}, dircont(j).name);
                        fnames{i} = fname;
                        % THIS ASSUMES SOME CONVENTION 
                        label{i} = dirj(5:(regexp(dirj,'MM')-2));
                        found  = true;
                    end
                end
            end

            if ~found
                warning(['no AVG or SUM file found in ' filepath{i}]);
            else
                info = imfinfo(fname);
                ims{i} = zeros([info(1).Height, info(1).Width 3], 'uint32');
            end

            % load data and determine new LUT
            %--------------------------------

            % only if no imposelim was passed
            warning('keeping all AVG images in memory, this could get out of control');

            % read RGB, not DIC, assuming DIC is first
            ci = 1;
            for c = uint16(corder{i}(2:4))

                imc = imread(fname, c);
                ims{i}(:,:,ci) = imc;

                imcmin = min(imc(:));
                imcmax = max(imc(:));

                % determine LUT limits per image
                slim = stretchlim(mat2gray(imc));
                minv = imcmin + (imcmax-imcmin)*slim(1);
                maxv = imcmin + (imcmax-imcmin)*slim(2);

                lims{i,ci} = double([minv maxv]);

                ci = ci+1;
            end 
        end
    end
    
    % combine LUT
    %----------------------------  
    if nargin == 4
        combinedLim = lims(1,:);

        for i = 2:numel(filepath)

            for c = 1:3 
                % as long as we use the same LUT, it should be fine
                % BUT: changing the lower limit, changes the relative values of
                % the background because it is like subtracting out the mean
                % and then looking at the ratio
                % its much nicer for us to look at in ilastik though
                %combinedLim{c}(1) = 0*min(combinedLim{c}(1), lims{i,c}(1));
                combinedLim{c}(1) = min(combinedLim{c}(1), lims{i,c}(1));
                combinedLim{c}(2) = max(combinedLim{c}(2), lims{i,c}(2));
            end
        end
    else
        disp('using external LUT');
        combinedLim = imposesLim;
    end
    
    % overwrite the lims so we use and return the combined one
    for i = 1:numel(filepath)
        for c = 1:3
            lims{i,c} = combinedLim{c};
        end
    end
    
    for j = 1:numel(batchIdx)
        
        for i = batchIdx{j}

            % save images with new LUT
            %----------------------------

            imScaled = zeros(size(ims{i}));
            for c= 1:3
                imScaled(:,:,c) = mat2gray(ims{i}(:,:,c), lims{i,c});
            end

            FD = imScaled;
            FD(:,:,3) = 0*FD(:,:,3);

            % save in each subdirectory
            if saveims
                imwrite(imScaled, fullfile(filepath{i}, [label{i} '_RGBscaled.tif']));    
                imwrite(imScaled(:,:,3), fullfile(filepath{i}, [label{i} '_nuclei.tif']));
                imwrite(FD, fullfile(filepath{i}, [label{i} '_FatDs.tif']));
                Ilim = lims(i,:);
                save(fullfile(filepath{i},'Ilim'), 'Ilim');
            end

            % save together in batch directory
            if saveims
                imwrite(imScaled(:,:,3), fullfile(datapath{j}, 'nuclei', [label{i} '_nuclei.tif']));
                imwrite(FD, fullfile(datapath{j}, 'FatDs', [label{i} '_FatDs.tif']));
            end

            % create subdirs for later results
            %----------------------------------

            if ~exist(fullfile(filepath{i}, segResultsDir), 'dir')
                mkdir(fullfile(filepath{i}, segResultsDir));
            end
            if ~exist(fullfile(filepath{i}, 'ilastik'), 'dir')
                mkdir(fullfile(filepath{i}, 'ilastik'));
            end
            if ~exist(fullfile(filepath{i}, 'analysis_results'), 'dir')
                mkdir(fullfile(filepath{i}, 'analysis_results'));
            end
        end
    end
    
    disp('preIlastik done');
end