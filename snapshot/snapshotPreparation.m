%--------------------------------------------------------------------------
%
%   produce images for Ilastik segmentation
%
%--------------------------------------------------------------------------

clear all;
close all;

% add subfolders of folder containing this script to path
%---------------------------------------------------------
[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

cd(testScriptPath);
cd('../..');
addpath(genpath(pwd));

% import parameters (data and code locations, ...)
%snapshotSegParametersLocal;
snapshotSegParametersLocal_exp020715;

%filepath = {'/Users/idse/Dropbox/Sprinzak/shared/snapshots 07.05.15/6h dox ilastic/dox 1h_5/'};
%corder = {[1 3 2 4]};

% check corder parameter
if numel(corder) ~= numel(filepath)
    warning('corder doesnt have enough elements, assuming default order of channels: 4 1 2 3');
    corder = {};
    for i = 1:numel(filepath), corder{i} = [4 1 2 3]; end
end

saveims = true;
[fnames, lims, flabel] = preIlastik(datapath, dataset, corder, saveims);
