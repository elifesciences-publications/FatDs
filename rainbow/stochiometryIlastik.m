%--------------------------------------------------------------------------
%
%   script for analyzing polarity
%
%--------------------------------------------------------------------------

clear all;
close all;

% add subfolders of folder containing this script to path
[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

cd(testScriptPath);
cd('..');
addpath(genpath(pwd));

datapath = '/Users/idse/Dropbox/Sprinzak/shared/images for fat-ds ratio calculations/analyzed 26.7.16';
cd(datapath);

files = {...
    fullfile(datapath, 'snaps fat ds 100x 5.5.16/44', '1_1000_W44.tif')...
    fullfile(datapath, 'snaps fat ds 100x 5.5.16/45', '1_1000_W45.tif')...
    fullfile(datapath, 'snaps fat ds 100x 5.5.16/47', '1_1000_W47.tif')...
    fullfile(datapath, 'snaps fat ds 100x 5.5.16/48', '1_1000_W48.tif')...
    fullfile(datapath, 'snaps fat ds 100x 5.5.16/49', '1_1000_W49.tif')...
    fullfile(datapath, 'snaps fat ds 100x 5.5.16/410', '1_1000_W410.tif')...
    fullfile(datapath, 'snaps fat ds 100x 5.5.16/411', '1_1000_W411.tif')...
    fullfile(datapath, 'snaps fat ds 100x 5.5.16/412', '1_1000_W412.tif')...
    fullfile(datapath, 'snaps fat ds 100x 5.5.16/415', '1_1000_W415.tif')...
    fullfile(datapath, 'snaps fat ds 100x 5.5.16/416', '1_1000_W416.tif')...
    fullfile(datapath, '24.7.16 100x co cul/418', '1_1000_W418.tif')...
    fullfile(datapath, '24.7.16 100x co cul/422', '1_1000_W422.tif')...
    fullfile(datapath, '24.7.16 100x co cul/423', '1_1000_W423.tif')...
    fullfile(datapath, '24.7.16 100x co cul/424', '1_1000_W424.tif')...
    fullfile(datapath, '24.7.16 100x co cul/425', '1_1000_W425.tif')...
    fullfile(datapath, '24.7.16 100x co cul/427', '1_1000_W427.tif')...
    fullfile(datapath, '24.7.16 100x co cul/428', '1_1000_W428.tif')...
    fullfile(datapath, '24.7.16 100x co cul/429', '1_1000_W429.tif')...
    fullfile(datapath, '24.7.16 100x co cul/430', '1_1000_W430.tif')...
    fullfile(datapath, '24.7.16 100x co cul/431', '1_1000_W431.tif')...
    fullfile(datapath, '24.7.16 100x co cul/432', '1_1000_W432.tif')...
    fullfile(datapath, '24.7.16 100x co cul/433', '1_1000_W433.tif')...
    fullfile(datapath, '24.7.16 100x co cul/434', '1_1000_W434.tif')...
    fullfile(datapath, '24.7.16 100x co cul/436', '1_1000_W436.tif')...
    fullfile(datapath, '24.7.16 100x co cul/437', '1_1000_W437.tif')...
    fullfile(datapath, '24.7.16 100x co cul/438', '1_1000_W438.tif')...
    fullfile(datapath, '24.7.16 100x co cul/439', '1_1000_W439.tif')...
    fullfile(datapath, '24.7.16 100x co cul/440', '1_1000_W440.tif')...
    fullfile(datapath, '24.7.16 100x co cul/441', '1_1000_W441.tif')...
    fullfile(datapath, '24.7.16 100x co cul/443', '1_1000_W443.tif')...
    fullfile(datapath, '24.7.16 100x co cul/444', '1_1000_W444.tif')...
    fullfile(datapath, '24.7.16 100x co cul/445', '1_1000_W445.tif')...
    fullfile(datapath, '24.7.16 100x co cul/448', '1_1000_W448.tif')...
    };

segmentationPath = {};
for i = 1:10
    segmentationPath{i} = '/Users/idse/Dropbox/Sprinzak/shared/images for fat-ds ratio calculations/analyzed 26.7.16/snaps fat ds 100x 5.5.16/ilastik segmentation';
end
for i = 11:numel(files)
    segmentationPath{i} = '/Users/idse/Dropbox/Sprinzak/shared/images for fat-ds ratio calculations/analyzed 26.7.16/24.7.16 100x co cul/ilastik segmentation';
end

zidx = {1,3,2,1,3,3,2,1,1,2,...
        3,2,2,1,2,2,2,2,1,2,2,1,1,3,1,3,2,1,1,1,2,3,3};

% % make a directory for saving the results
% resultsdir = fullfile(filepath, [barefname '_analysis']);
% if ~exist(resultsdir, 'dir')
%     mkdir(resultsdir);
% end

% SETTINGS! 
%-----------------

% channel indices
DsC = 1;
FatC = 2;

% % resolution, nm per pixel
% res = 60;

% SETTINGS
thickness = 11;

% emissions peak wavelengths, 530, 610

%% print out resolution of each file

allRes = {};

for fi = 1:numel(files)

    %-----------------------------
    % read the data
    %-----------------------------

    [filepath,barefname,ext] = fileparts(files{fi});
    resultsdir = filepath;
    fname = [barefname ext]; 
    [data meta] = readStack(filepath, fname);
    data = squeeze(data);
    allRes{fi} = round(1000*meta.xres); % resolution in nm/pixel
end

diary(fullfile(datapath,'resolutions.txt'));
for fi = 1:numel(files)
    disp(['file number: ' num2str(fi)]);
    disp(['file name: ' files{fi}]);
    disp(['resolution: ' num2str(allRes{fi})]);
    disp('-');
end
diary off

%% get boundary intensities 

bdryFatTot = [];
bdryDsTot = [];
fatBdry = {};
dsBdry = {};

bdryDs = zeros([1 numel(files)]);
bdryFat = zeros([1 numel(files)]);

for fi = 1:numel(files)

    %-----------------------------
    % read the data
    %-----------------------------

    [filepath,barefname,ext] = fileparts(files{fi});
    resultsdir = filepath;
    fname = [barefname ext]; 
    [data meta] = readStack(filepath, fname);
    data = squeeze(data);
%     if fi < 7
%         data = squeeze(data(:,:,zidx{fi},[3 2 1]));
%     end

    res = round(1000*meta.xres); % resolution in nm/pixel

    % read Ilastik segmentation
    segfname = dir([fullfile(segmentationPath{fi}, barefname) '*']);
    segfname = segfname.name;
    bdryMask = imread(fullfile(segmentationPath{fi}, segfname)) == 2;
    %bdryMask = imclose(bdryMask,strel('disk',2));
    
    Ds = data(:,:,DsC);
    Fat = data(:,:,FatC);
    
    Ds = Ds - imopen(Ds,strel('disk',100));
    Fat = Fat - imopen(Fat,strel('disk',100));
    
    bdryDsMean(fi) = mean(Ds(bdryMask));
    bdryFatMean(fi) = mean(Fat(bdryMask)); 

    bdryDsSum(fi) = sum(Ds(bdryMask));
    bdryFatSum(fi) = sum(Fat(bdryMask)); 
%     % overlay image
%     bdryEdge = bdryMask - imerode(bdryMask, strel('disk',1));
%     R = imadjust(Ds);% + 255*uint8(DsEdge);
%     G = imadjust(Fat);% + 255*FatEdge;
%     R = mat2gray(Ds);% + 255*uint8(DsEdge);
%     G = mat2gray(Fat);% + 255*FatEdge;
%     B = (2^16-1)*bdryEdge;
%     R(bdryEdge > 0) = 0;
%     G(bdryEdge > 0) = 0;
%     overlay = cat(3, R, G, B);
%     
%     fname = [barefname '_segoverlay_unsaturated.tif'];
%     imwrite(overlay, fullfile(segmentationPath,fname));
%     
    %imshow(overlay)
end

%% plot bdry Fat vs Ds

figure, 
plot(bdryDsSum,bdryFatSum,'x')
% hold on 
% for i = 1:numel(bdryDsSum)
%     text(bdryDsSum(i),bdryFatSum(i),num2str(i))
% end
% hold off
xlabel('sum bdry Ds');
ylabel('sum bdry Fat');
C = corrcoef(bdryDsSum,bdryFatSum);
title(['correlation: ' num2str(C(2),3)])
saveas(gcf,fullfile(datapath,'bdryFatDsSum.png'));
saveas(gcf,fullfile(datapath,'bdryFatDsSum.fig'));
% MEAN VS SUM?

%% plot bdry Fat vs Ds

figure, 
plot(bdryDsMean,bdryFatMean,'x')
% hold on 
% for i = 1:numel(bdryDsMean)
%     text(bdryDsMean(i),bdryFatMean(i),num2str(i))
% end
% hold off
xlabel('mean bdry Ds');
ylabel('mean bdry Fat');
C = corrcoef(bdryDsMean, bdryFatMean);
title(['correlation: ' num2str(C(2),3)])
saveas(gcf,fullfile(datapath,'bdryFatDsMean.png'));
saveas(gcf,fullfile(datapath,'bdryFatDsMean.fig'));
% MEAN VS SUM?
