%---------------------------------------------------------
% SCRIPT TO MAKE AVERAGE OF NUMBER OF SLICES OF Z-STACK
%---------------------------------------------------------

clear all;
close all;

% add subfolders of folder containing this script to path
%---------------------------------------------------------
[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

cd(testScriptPath);
cd('../..');
addpath(genpath(pwd));

% OLGA! parameters
%---------------------------------------------------------

snapshotSegParametersLocal;

zi = 3; % first slice to include in average
zf = 7; % last slice to include in average

%regenerate filepath and batchIdx
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

%scan directories voor data
for j = 1:numel(batchIdx)

    for i = batchIdx{j}

        found = false;

        %identify file starting with SUM
        %----------------------------

        dircont = dir(filepath{i});
        for j = 1:numel(dircont)

            dirj = dircont(j).name;

            if ~isempty(regexp(dirj,'ome.tif','ONCE')) && ...
                isempty(regexp(dirj,'AVG','ONCE'))
            
                fnames{i} = dircont(j).name;
                found  = true;
                
                makeAVG(filepath{i},fnames{i},zi,zf);
            end
        end
        if ~found 
            warning(['no ome.tif file found in ' filepath{i}]);
        end
    end
end