%--------------------------------------------------------------------------
%
%   script for analyzing movies
%
%--------------------------------------------------------------------------

clear all;
close all;

% add subfolders of folder containing this script to path
[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

cd(testScriptPath);
cd('..');
addpath(genpath(pwd));

datapath = {};
barefnameFormat = {};
tmin = {};
tmax = {};

datapath{1} = '/Users/idse/Dropbox/Sprinzak/shared/movies/analysis 1.8.13 movie 3';
barefnameFormat{1} = 'z8-10-%.4d';
tmin{1} = 2;
tmax{1} = 19;

datapath{2} = '/Users/idse/Dropbox/Sprinzak/shared/movies/18.11.13 movie 2';
barefnameFormat{2} = 'AVG_002_z3-1-%.4d';
tmin{2} = 10;
tmax{2} = 100;

datapath{3} = '/Users/idse/Dropbox/Sprinzak/shared/movies/4.3.14/4';
barefnameFormat{3} = 'AVG_Composite4-1-%.4d';
tmin{3} = 11;
tmax{3} = 77;

datapath{4} = '/Users/idse/Dropbox/Sprinzak/shared/movies/4.3.14/18';
barefnameFormat{4} = 'AVG_Composite18-RGB-%.4d';
tmin{4} = 5;
tmax{4} = 24;

%% extract data

movieIdx = 2;

outputDir = fullfile(datapath{movieIdx}, 'analysis_results');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

fnameFormat = [barefnameFormat{movieIdx} '.tif'];
segFnameFormat = [barefnameFormat{movieIdx} '_seg.tif'];
bdrysegFnameFormat = [barefnameFormat{movieIdx} '_bdryseg.tif'];
T = tmax{movieIdx};

% levels
TotFatLevel = zeros([1 T]);
TotDsLevel = zeros([1 T]);
bdryFatLevel = zeros([1 T]);
bdryDsLevel = zeros([1 T]);
cellFatLevel = zeros([1 T]);
cellDsLevel = zeros([1 T]);

% background
Dsbg = zeros([1 T]);
Fatbg = zeros([1 T]);

% areas
TotDsA = zeros([1 T]);
cellDsA = zeros([1 T]);
TotFatA = zeros([1 T]);
cellFatA = zeros([1 T]);
bdryA = zeros([1 T]);

%%

FMIP = zeros([512 512 T]);%,'uint16');
DMIP = zeros([512 512 T]);%,'uint16');
for i = 1:T
    FMIP(:,:,i) = imread(fullfile(datapath{movieIdx},'C2-SUM_002.tif'),i);
    DMIP(:,:,i) = imread(fullfile(datapath{movieIdx},'C1-SUM_002.tif'),i);
end

%%

for i = 1:T
    
    %data = imread(fullfile(datapath{movieIdx},sprintf(fnameFormat,i)));
    
    seg = imread(fullfile(datapath{movieIdx},sprintf(segFnameFormat,i)));
    bdryseg = imread(fullfile(datapath{movieIdx},sprintf(bdrysegFnameFormat,i)));
    
    %Fat = data(:,:,2);
    %Ds = data(:,:,1);
    Fat = FMIP(:,:,i);
    Ds = DMIP(:,:,i);
    
    % clean up teh segmentation a little
    FatMask = imclose(seg == 2, strel('disk',1));
    DsMask = imclose(seg == 1, strel('disk',1));
    bdryMask = bdryseg == 1;
    totalMask = FatMask | DsMask | bdryMask;
    bgMask = imdilate(totalMask, strel('disk',30)) - totalMask > 0;
    
    FatMask = bwareaopen(FatMask, 100);
    DsMask = bwareaopen(DsMask, 100);
    
    % visualize the segmentation and save overlay
    bdryEdge = bdryMask - imerode(bdryMask, strel('disk',2));
    FatEdge = FatMask - imerode(FatMask, strel('disk',2));
    DsEdge = DsMask - imerode(DsMask, strel('disk',2));
    R = Ds + 255*DsEdge;
    G = double(Fat) + 255*FatEdge;
    B = 255*bdryEdge;
    R(bdryEdge > 0) = 0;
    G(bdryEdge > 0) = 0;
    overlay = cat(3, R, G, B);
    %imshow(overlay,[])

    [~,x,~] = fileparts(sprintf(fnameFormat,i));
    fname = [x '_segoverlay.tif'];
    imwrite(overlay, fullfile(outputDir,fname));
    
    % background 
    Dsbg(i) = median(Ds(bgMask));
    Fatbg(i) = median(Fat(bgMask));
    
    % areas
    TotFatA(i) = sum(sum(FatMask | bdryMask));
    TotDsA(i) = sum(sum(DsMask | bdryMask));
    cellDsA(i) = sum(sum(DsMask & ~bdryMask));
    cellFatA(i) = sum(sum(FatMask & ~bdryMask));
    bdryA(i) = sum(bdryMask(:));
    
    % read out the levels
    bdryFatLevel(i) = sum(Fat(bdryMask));
    bdryDsLevel(i) = sum(Ds(bdryMask));
    TotFatLevel(i) = sum(Fat(FatMask | bdryMask));
    TotDsLevel(i) = sum(Ds(DsMask| bdryMask));
    cellFatLevel(i) = sum(Fat(FatMask & ~bdryMask));
    cellDsLevel(i) = sum(Ds(DsMask & ~bdryMask));
    
    %bgFatLevel(i) = mean(Fat(~(FatMask | DsMask | bdryMask)));
    %bgDsLevel(i) = mean(Ds(~(FatMask | DsMask | bdryMask)));
end


%%

tidx = 1:40;%tmin{movieIdx}:tmax{movieIdx};

TotFatbgsub = TotFatLevel - TotFatA.*Fatbg;
TotDsbgsub = TotDsLevel - TotDsA.*Dsbg;
bdryFatbgsub = bdryFatLevel - bdryA.*Fatbg;
bdryDsbgsub = bdryDsLevel - bdryA.*Dsbg;

% bdryFatbgsub(bdryFatbgsub <= 0)=0.1;
bdryDsbgsub(bdryDsbgsub <= 0)=0.1;

TotFatbgsub = TotFatbgsub(tidx);
TotDsbgsub = TotDsbgsub(tidx);
bdryFatbgsub = bdryFatbgsub(tidx);
bdryDsbgsub = bdryDsbgsub(tidx);

ext = {'.png','.fig'};
axislim = [tidx(1) tidx(end) 10^3 10^8];

for i = 1%:numel(ext)
    
    figure,
    semilogy(tidx, TotFatbgsub)
    axis(axislim)
    %saveas(gcf, fullfile(outputDir,['TotFatLevel' ext{i}]));
    
    figure,
    semilogy(tidx,TotDsbgsub)
    axis(axislim)
    %saveas(gcf, fullfile(outputDir,['TotDsLevel' ext{i}]));

%     plot(cellFatLevel)
%     saveas(gcf, fullfile(outputDir,['cellFatLevel' ext{i}]));
%     
%     plot(cellDsLevel)
%     saveas(gcf, fullfile(outputDir,['cellDsLevel' ext{i}]));
%     
    figure,
    semilogy(tidx,bdryFatbgsub)
    axis(axislim)
    %saveas(gcf, fullfile(outputDir,['bdryFatLevel' ext{i}]));
    
    figure,
    semilogy(tidx,bdryDsbgsub)
    axis(axislim)
    %saveas(gcf, fullfile(outputDir,['bdryDsLevel' ext{i}]));
    
    
    loglog(bdryFatbgsub,bdryDsbgsub,'x')
    %saveas(gcf, fullfile(outputDir,['bdryFatVsDs' ext{i}]));
    
%     
%     plot(TotDsLevel, bdryDsLevel)
%     saveas(gcf, fullfile(outputDir,['bdryVsTotDsLevel' ext{i}]));
% 
%     loglog(TotDsLevel, bdryDsLevel)
%     saveas(gcf, fullfile(outputDir,['bdryVsTotDsLevelLog' ext{i}]));
end

%save(fullfile(outputDir, 'bdryIntensities'),'bdryFatbgsub','bdryDsbgsub')

%% combined plot stochiometry of movies and fixed

colors = lines(4);
Nall = 10^2;
figure('Position',[1 1 800 800])
name = {};

for movieIdx = 1:4
    
    outputDir = fullfile(datapath{movieIdx}, 'analysis_results');
    tidx = tmin{movieIdx}:tmax{movieIdx};
    m = load(fullfile(outputDir, 'bdryIntensities'));
    
    bdryDsNormalized = m.bdryDsbgsub*(mean(m.bdryFatbgsub)/mean(m.bdryDsbgsub));
    %bdryDsNormalized = bdryDsNormalized./mean(bdryDsNormalized);
    
    bdryFatNormalized = m.bdryFatbgsub;
    %bdryDsNormalized = m.bdryDsbgsub./m.bdryDsbgsub(1);
    
    bdryDsNormalized = bdryDsNormalized./Nall;
    bdryFatNormalized = bdryFatNormalized./Nall;
    
    if movieIdx > 1
        hold on
    end
    loglog(bdryFatNormalized, bdryDsNormalized,'-x','Color',colors(movieIdx,:))
    %loglog(m.bdryFatbgsub, m.bdryDsbgsub,'xr')
    axis([0.5 10^4 0.5 10^4])
    
    bla = strsplit(datapath{movieIdx},'/');
    name{movieIdx} = bla{end};
    %clf
    %plot(log(bdryFatNormalized), log(bdryDsNormalized),'-xr')
end
legend(name,'Location','NorthWest');
hold off

saveas(gcf, '/Users/idse/Dropbox/Sprinzak/shared/movies/moviesCombinedPlot.png');

