function [fnames, lims, label] = postIlastik(datapath, dataset)
    % move segmentations back to the right place

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
                error(['no AVG or SUM file found in ' filepath{i}]);
            end
        end
    end
    
    for j = 1:numel(batchIdx)
        
        for i = batchIdx{j}
            
            S = load(fullfile(filepath{i},'Ilim'));
            for c = 1:3
                lims{i,c} = S.Ilim{c};
            end
            
            f1 = fullfile(datapath{j}, 'nuclei', [label{i} '_nuclei_seg.tif']);
            f2 = fullfile(filepath{i}, 'ilastik', [label{i} '_nuclei_seg.tif']);
            
            if exist(f1,'file') %&& ~exist(f2,'file')
                copyfile(f1, f2);
            end
            
            f1 = fullfile(datapath{j}, 'FatDs', [label{i} '_FatDs_seg.tif']);
            f2 = fullfile(filepath{i}, 'ilastik', [label{i} '_FatDs_seg.tif']);
            
            if exist(f1,'file') %&& ~exist(f2,'file')
                copyfile(f1, f2);
            end
            
            f1 = fullfile(datapath{j}, 'FatDs', [label{i} '_FatDs_bdryseg.tif']);
            f2 = fullfile(filepath{i}, 'ilastik', [label{i} '_FatDs_bdryseg.tif']);
            
            if exist(f1,'file') %&& ~exist(f2,'file')
                copyfile(f1, f2);
            end
        end
    end
    
    disp('postIlastik done');
end