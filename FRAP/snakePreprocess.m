%--------------------------------------------------------------------------
%
%   script for preprocessing FRAP data
%
%--------------------------------------------------------------------------

clear all;
close all;

% add subfolders of folder containing this script to path
[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(testScriptPath));

datapath = '/Users/idse/David_visit_data/';
cd(datapath);
[fname, filepath] = uigetfile('*.tif');

[~,barefname,~] = fileparts(fname);

% channel indices
DsC = 2;
FatC = 3;

%%
%-----------------------------
% read the data
%-----------------------------

[data nChannels] = readStack(filepath, fname);
nTimePts = size(data,5);


%% 
%-----------------------------
% enhance for active contours
%-----------------------------

zSlice = 1;
forBdryDetect = zeros([size(data,1), size(data,2), nTimePts]);

for t = 1:nTimePts
    
    slice = squeeze(data(:,:,zSlice, FatC, t));
    th = imtophat(slice, strel('square',4)); 
    forBdryDetect(:,:,t) = mat2gray(th);
end

writeTiff(forBdryDetect, fullfile(filepath, ['forBdryDetect' barefname]));

%% try what works

zSlice = 1;
t = 2;
slice = squeeze(data(:,:,zSlice, FatC, t));
slice = slice - imopen(slice, strel('square', 4));
slice = mat2gray(slice);

%% find boundaries without snake

thresh = mean(slice(:)) + 2*std(slice(:));
imshow(bwareaopen(slice > thresh, 100))