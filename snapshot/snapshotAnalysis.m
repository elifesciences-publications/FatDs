%--------------------------------------------------------------------------
%
%   analysis of snapshot segmentation
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

% import parameters:
% defines filepath cell array, batchLabels cell array and saveIntermediates
%snapshotSegParametersLocal_exp280715;
snapshotSegParametersLocal_exp020715
%snapshotSegParametersLocal_hysteresis

% check corder parameter
if numel(corder) ~= numel(filepath)
    warning('corder doesnt have enough elements, assuming default order of channels: 4 1 2 3');
    corder = {};
    for i = 1:numel(filepath), corder{i} = [4 1 2 3]; end
end

[fnames, lims, flabel] = postIlastik(datapath, dataset);

% read cellLayer objects and other segmentation results
CL = {};
bdryI = {};
nImages = numel(filepath);
nBatches = numel(batchLabel);

% for robustness of batches
goodseg = [];

for i = 1:nImages
    
    % PAY ATTENTION: there was _new at the end here at some point
    % _seg_exclNuc = exclude accumulating boundaries (tired in naming)
    fname = fullfile(filepath{i},'matlab_seg',[flabel{i} '_seg_new']);
    if exist([fname '.mat'],'file')
        S = load(fname);
        CL{i} = S.cellLayer;
%         bdryI{i} = struct('bdryFI', {S.bdryFI}, 'bdryDI', {S.bdryDI},...
%                             'asBdry', {S.assignedBdry});
        goodseg = [goodseg i];
    else
        warning([fname ' doesnt exist']);
    end
end

for i = 1:nBatches
    batchIdx{i} = intersect(batchIdx{i}, goodseg);
end

% save results
%--------------------------------------------------
saveResults = false;

% PAY ATTENTION: there was _new at the end here at some point
ext = '526.png';

stats = {};

% fixed samples, not time
t = 1;

Nbins = 50;

%%
%--------------------------------------------------------------
% cell intensity distributions
%--------------------------------------------------------------

for i = goodseg
    
    cellLayer = CL{i};
    
    % cells
    cellStates = cat(1,cellLayer.cells{t}.state);

    FatIdx = cellStates(:,5)==1;
    DsIdx = cellStates(:,6)==1;
    otherIdx = cellStates(:,5)~=1 & cellStates(:,6)~=1;

    % intensities
    IFat = cellStates(:,1);
    IDs = cellStates(:,2);
    mean(IDs)

    % bins
    Ibins = {};
    adjust = [0.7 1 1];
    for c = 1:3
        Imax = lims{i,c}(2)*adjust(c);
        Imin = lims{i,c}(1);
        Irange = Imax - Imin;
        Ibins{c} = Imin:Irange/Nbins:Imax;
    end
    
    % index of cells that have a boundary
    cellsWBdry = false(size(FatIdx));
    % index of cells that have no boundary but do have a Fat-Ds interface
    cellsWFDif = false(size(FatIdx)); 
    
    for ci = 1:numel(cellsWBdry)
        
        % bond states in cell
        cbi = cellLayer.cells{t}(ci).bondInd;
        bsOfCell = cat(1, cellLayer.bonds{t}(cbi).state);
        
        % 4th column of bond.state is boolean bdry
        if ~isempty(bsOfCell) && any(bsOfCell(:,4)==1)
            cellsWBdry(ci)=true;
        end 
        
        if ~isempty(bsOfCell) && any(bsOfCell(:,1)==1)
            cellsWFDif(ci)=true;
        end 
    end

    % FDIF == SEG & FDIF
    % cells that have no boundary but do have a Fat-Ds interface
    % to check against cells that have no bdry in general -> same
    stats{i}.IDsFDIF = makeDist(IDs(cellsWFDif & DsIdx), Ibins{1});
    stats{i}.IFatFDIF = makeDist(IFat(cellsWFDif & FatIdx), Ibins{2}); 

    % BDRY == SEG & BDRY
    stats{i}.IDsBdry = makeDist(IDs(cellsWBdry & DsIdx), Ibins{1});
    stats{i}.IFatBdry = makeDist(IFat(cellsWBdry & FatIdx), Ibins{2}); 

    % NO BDRY == SEG & NO BDRY
    stats{i}.IDsNoBdry = makeDist(IDs(~cellsWBdry & DsIdx), Ibins{1});
    stats{i}.IFatNoBdry = makeDist(IFat(~cellsWBdry & FatIdx), Ibins{2}); 

    % NEW: logscale binning
    
%     logbinsDs = exp(linspace(log(min(Ibins{1})), log(max(Ibins{1})),50));
%     stats{i}.IDsBdryLogBin = makeDist(IDs(cellsWBdry & DsIdx)-100, logbinsDs);
%     stats{i}.IDsNoBdryLogBin = makeDist(IDs(~cellsWBdry & DsIdx)-100, logbinsDs);

%     logbinsFat = exp(linspace(log(min(Ibins{2})), log(max(Ibins{2})), Nbins));
%     stats{i}.IFatBdryLogBin = makeDist(IFat(cellsWBdry & FatIdx), logbinsFat);
%     stats{i}.IFatNoBdryLogBin = makeDist(IFat(~cellsWBdry & FatIdx), logbinsFat);
  
    % stats after substracting background value, lambdaI cutoff 
    lambdaI = 10;
    IDsmin = Ibins{1}(1);
    IDsmax = Ibins{1}(end);
    logbinsDs = exp(linspace(log(lambdaI), log(IDsmax - IDsmin), Nbins));
    stats{i}.IDsBdryLogBin = makeDist(IDs(cellsWBdry & DsIdx) - IDsmin, logbinsDs);
    stats{i}.IDsNoBdryLogBin = makeDist(IDs(~cellsWBdry & DsIdx)  - IDsmin, logbinsDs);
    
    IFatmin = Ibins{2}(1);
    IFatmax = Ibins{2}(end);
    logbinsFat = exp(linspace(log(lambdaI), log(IFatmax - IFatmin), Nbins));
    stats{i}.IFatBdryLogBin = makeDist(IFat(cellsWBdry & FatIdx) - IFatmin, logbinsFat);
    stats{i}.IFatNoBdryLogBin = makeDist(IFat(~cellsWBdry & FatIdx)  - IFatmin, logbinsFat);
    
    % ALL: intensity distribution of all cells
    stats{i}.IDs = makeDist(IDs, Ibins{1});
    stats{i}.IFat = makeDist(IFat, Ibins{2});

    % SEG: intensity distribution of segmented cells
    stats{i}.IDsSeg = makeDist(IDs(DsIdx), Ibins{1});
    stats{i}.IFatSeg = makeDist(IFat(FatIdx), Ibins{2});
end

%%
%--------------------------------------------------------------
% combine cell intensity of each batch
%--------------------------------------------------------------

statsTot = {};

for i = 1:nBatches
    
    IDsTot = [];
    IFatTot = [];
    
    IDsFDIFTot = [];
    IFatFDIFTot = [];
    
    IDsBdryTot = [];
    IFatBdryTot = [];
    IDsNoBdryTot = [];
    IFatNoBdryTot = [];
    IDsSegTot = [];
    IFatSegTot = [];
    batchStat = struct();
    
    % combine data
    %-----------------------------
    for j = 1:numel(batchIdx{i})
        
          idx = batchIdx{i}(j);
          
          % ALL
          IDsTot = cat(1, IDsTot, stats{idx}.IDs.table);
          IFatTot = cat(1, IFatTot, stats{idx}.IFat.table);

          % SEG FDIF
          IDsFDIFTot = cat(1, IDsFDIFTot, stats{idx}.IDsFDIF.table);
          IFatFDIFTot = cat(1, IFatFDIFTot, stats{idx}.IFatFDIF.table);

          % SEG BDRY
          IDsBdryTot = cat(1, IDsBdryTot, stats{idx}.IDsBdry.table);
          IFatBdryTot = cat(1, IFatBdryTot, stats{idx}.IFatBdry.table);

          % SEG NO BDRY
          IDsNoBdryTot = cat(1, IDsNoBdryTot, stats{idx}.IDsNoBdry.table);
          IFatNoBdryTot = cat(1, IFatNoBdryTot, stats{idx}.IFatNoBdry.table);

          % SEG
          IDsSegTot = cat(1, IDsSegTot, stats{idx}.IDsSeg.table);
          IFatSegTot = cat(1, IFatSegTot, stats{idx}.IFatSeg.table);
    end
    
    % make combined distribution
    %-----------------------------
    
    % ALL
    batchStat.IDs = makeDist(IDsTot, Ibins{1});
    batchStat.IFat = makeDist(IFatTot, Ibins{2});
     
    % SEG
    batchStat.IDsSeg = makeDist(IDsSegTot, Ibins{1});
    batchStat.IFatSeg = makeDist(IFatSegTot, Ibins{2});
    
    % SEG FDIF LOG
    batchStat.IDsFDIFLogBin = makeDist(IDsFDIFTot, logbinsDs);
    batchStat.IFatFDIFLogBin = makeDist(IFatFDIFTot, logbinsFat);
    
    % SEG NO BDRY
    batchStat.IDsNoBdry = makeDist(IDsNoBdryTot, Ibins{1});
    batchStat.IFatNoBdry = makeDist(IFatNoBdryTot, Ibins{2});

    % SEG NO BDRY LOG
    batchStat.IDsNoBdryLogBin = makeDist(IDsNoBdryTot - IDsmin, logbinsDs);
    batchStat.IFatNoBdryLogBin = makeDist(IFatNoBdryTot - IDsmin, logbinsFat);
          
    % SEG BDRY
    batchStat.IDsBdry = makeDist(IDsBdryTot, Ibins{1});
    batchStat.IFatBdry = makeDist(IFatBdryTot, Ibins{2});

    % SEG BDRY LOG
    batchStat.IDsBdryLogBin = makeDist(IDsBdryTot  - IDsmin, logbinsDs);
    batchStat.IFatBdryLogBin = makeDist(IFatBdryTot - IDsmin, logbinsFat);
          
    statsTot{i} = batchStat;
end

%% check that distributions are made correctly and compare 2D dist

bi = 12;

binstmp = 100:100:2000;
dist = hist(IDsNoBdryTot - 100, binstmp);
dist = dist./sum(dist(:));
% 
% figure, plot(binstmp, dist,'-x')

figure,
%semilogx(logbinsDs(1:end-1), statsTot{bi}.IDsFDIFLogBin.dist(1:end-1), '-x', 'Color', 'b', 'LineWidth', 2);
% as expected cells with an interface to the other kind (~2/3) 
% have the same distribution of intensities
semilogx(logbinsDs(1:end-1), statsTot{bi}.IDsNoBdryLogBin.dist(1:end-1), '-x', 'Color', 'r', 'LineWidth', 2);
hold on
semilogx(logbinsDs(1:end-1), statsTot{bi}.IDsBdryLogBin.dist(1:end-1), '-x', 'Color', 'b', 'LineWidth', 2);
%semilogx(logbinsDs(1:end-1), log(logbinsDs(1:end-1))/300,'-x');
%semilogx(logbinsDs(1:end-1), dist(1:end-1), '-x', 'Color', 'g', 'LineWidth', 2);
xlim([logbinsDs(1) logbinsDs(end)]);
hold off

figure,
plot(Ibins{1}(1:end-1), statsTot{bi}.IDsNoBdry.dist(1:end-1), '-x', 'Color', 'r', 'LineWidth', 2);
hold on
plot(Ibins{1}(1:end-1), statsTot{bi}.IDsBdry.dist(1:end-1), '-x', 'Color', 'b', 'LineWidth', 2);
xlim([Ibins{1}(1) Ibins{1}(end)]);
hold off

figure,
semilogx(Ibins{1}(1:end-1), statsTot{bi}.IDsNoBdry.dist(1:end-1), '-x', 'Color', 'r', 'LineWidth', 2);
hold on
semilogx(Ibins{1}(1:end-1), statsTot{bi}.IDsBdry.dist(1:end-1), '-x', 'Color', 'b', 'LineWidth', 2);
xlim([Ibins{1}(1) Ibins{1}(end)]);
hold off


%% look at Fat avg over time and Fat median

M = [];
for i = 1
    for j = 1:numel(batchIdx{i})
        idx = batchIdx{i}(j);
        M = cat(1, M, stats{idx}.IFatSeg.table);
        [j mean(stats{idx}.IFatSeg.table)]
    end
end

FatAvg = [];
FatMedian = [];
for i = 1:nBatches
    
    FatAvg(i) = mean(statsTot{i}.IFatSeg.table);
    FatMedian(i) = median(statsTot{i}.IFatSeg.table);
    
    %FatAvg(i) = mean(statsTot{i}.IFatBdry.table);
    %FatMedian(i) = median(statsTot{i}.IFatBdry.table);
end

figure,
plot(batchTimes, FatMedian,'-x','LineWidth',2)
%ylim([180 240])
hold on
plot(batchTimes, FatAvg,'-x','LineWidth',2)
hold off
legend({'median','mean'});
title('Fat intensity in Fat cells')
%title('Fat intensity in Fat cells with boundary')

saveResults = true;
if saveResults
    
    saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['FatVsTime' ext]));    
end

%%
% visualize Ds 

ncols = 2;

scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)/2 scrsz(3)*0.7 scrsz(4)]);

plotTitles = {'Ds I,', 'Ds I in Ds cells,',...
                'Ds I in cells w/ bdry,', 'Ds I in cells w/o bdry,'};
fieldnames = {'IDs','IDsSeg','IDsBdry', 'IDsNoBdry'};
nrows = numel(plotTitles);

for row = 1:nrows
    plotInfo = struct('plotTitle', plotTitles{row});
    subplotDist(nrows, ncols, row, plotInfo, stats, fieldnames{row}, flabel,...
                    batchIdx, batchLabel)
    %subplotDist(nrows, ncsols, row, plotInfo, stats(batchIdx{3}), fieldnames{row}, flabel(batchIdx{3}));
end
set(gcf, 'Color', 'w');
if saveResults
    
    saveas(h, fullfile(combinedOutputPath, 'analysis_results', ['DsIntensityAll' ext]));
    
    % also compare individual images within batch
    for i = 1:nBatches
        clf 
        if ~isempty(batchIdx{i})
            for row = 1:nrows
                plotInfo = struct('plotTitle', plotTitles{row});
                subplotDist(nrows, ncols, row, plotInfo, stats(batchIdx{i}),...
                                fieldnames{row}, flabel(batchIdx{i}));
            end
        end
        saveas(h, fullfile(datapath{i}, ['DsIntensityAll_' batchLabel{i} ext]));
    end
    
    % and the combined distributions
    clf
    for row = 1:nrows
        plotInfo = struct('plotTitle', plotTitles{row});
        subplotDist(nrows, ncols, row, plotInfo, statsTot, fieldnames{row}, batchLabel)
    end
    saveas(h, fullfile(combinedOutputPath, 'analysis_results', ['DsIntensityCombined' ext]));
     close
end

%% only 20h average for figure

figure,
batchLabel{bi}

ext = '.fig';
clf
for row = 1:nrows
    plotInfo = struct('plotTitle', plotTitles{row});
    subplotDistSingleTimeForFigure(nrows, ncols, row, plotInfo, statsTot(bi), fieldnames{row}, batchLabel(bi));
end
%saveas(h, fullfile(combinedOutputPath, 'analysis_results', ['DsIntensity20hAvg' ext]));

%% only 20 h combined distribution Ds for figure

clf
%h = gcf;
plot(statsTot{bi}.IDsNoBdry.bins, statsTot{bi}.IDsNoBdry.dist, 'Color', 'r', 'LineWidth', 2);
hold on
plot(statsTot{bi}.IDsBdry.bins, statsTot{bi}.IDsBdry.dist, 'Color', 'k', 'LineWidth', 2);
hold off
legend({'Ds no boundary','Ds boundary'})
%saveas(h, fullfile(combinedOutputPath, 'analysis_results', ['DsIntensity20hAvgCombined' ext]));

%% only 20 h combined distribution Ds on log scale for figure

figure,
clf
h = gcf;

% bin def:
% logbinsDs = exp(linspace(log(min(Ibins{1})), log(max(Ibins{1})),50));

bins = statsTot{bi}.IDsNoBdryLogBin.bins(1:end-1);
nobdrydist = statsTot{bi}.IDsNoBdryLogBin.dist(1:end-1);
bdrydist = statsTot{bi}.IDsBdryLogBin.dist(1:end-1);

[x,y] = histForBarlikePlot(bins, nobdrydist');
semilogx(x, y, '-', 'Color', 'r', 'LineWidth', 2);
hold on
[x,y] = histForBarlikePlot(bins, bdrydist');
semilogx(x, y, '-', 'Color', 'k', 'LineWidth', 2);
hold off
legend({'Ds no boundary','Ds boundary'})
axis([10 bins(end) 0 0.05])
%axis([10^2 1.2*10^3 0 0.05])
if saveResults
    saveas(h, fullfile(combinedOutputPath, 'analysis_results', ['DsIntensity20hAvgCombinedLogBinnedX' ext]));
end

%% only 20 h combined distribution Fat for figure

clf
plot(statsTot{bi}.IFatNoBdry.bins, statsTot{bi}.IFatNoBdry.dist, 'Color', 'r', 'LineWidth', 2);
hold on
plot(statsTot{bi}.IFatBdry.bins, statsTot{bi}.IFatBdry.dist, 'Color', 'k', 'LineWidth', 2);
hold off
legend({'Fat no boundary','Fat boundary'})
%saveas(h, fullfile(combinedOutputPath, 'analysis_results', ['FatIntensity20hAvgCombined' ext]));

%% only 20 h combined distribution Fat on log scale for figure

figure

bins = statsTot{bi}.IFatNoBdryLogBin.bins(1:end-1);
nobdrydist = statsTot{12}.IFatNoBdryLogBin.dist(1:end-1);
bdrydist = statsTot{12}.IFatBdryLogBin.dist(1:end-1);
[x,y] = histForBarlikePlot(bins, nobdrydist');
semilogx(x, y, '-', 'Color', 'r', 'LineWidth', 2);
hold on
[x,y] = histForBarlikePlot(bins, bdrydist');
semilogx(x, y, '-', 'Color', 'k', 'LineWidth', 2);
hold off
axis([10 bins(end) 0 0.1])
legend({'Fat no boundary','Fat boundary'})
if saveResults
    saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['FatIntensity20hAvgCombinedLogBinned' ext]));
end

%%
% visualize Fat 

scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)/2 scrsz(3)*0.7 scrsz(4)]);

plotTitles = {'Fat intensity', 'Fat intensity in cells labeled as Fat',...
    'Fat intensity in cells with bdry,', 'Fat I w/o bdry,'};
fieldnames = {'IFat','IFatSeg','IFatBdry','IFatNoBdry'};
nrows = numel(plotTitles);

for row = 1:nrows
    plotInfo = struct('plotTitle', plotTitles{row});
    subplotDist(nrows, ncols, row, plotInfo, stats, fieldnames{row}, flabel,...
                    batchIdx, batchLabel)
end

set(gcf, 'Color', 'w');
if saveResults
    
    saveas(h, fullfile(combinedOutputPath, 'analysis_results', ['FatIntensityAll' ext]));
    
    % also compare individual images within batch
    for i = 1:nBatches
        clf 
        for row = 1:nrows
            plotInfo = struct('plotTitle', plotTitles{row});
            subplotDist(nrows, ncols, row, plotInfo, stats(batchIdx{i}), fieldnames{row}, flabel(batchIdx{i}));
        end
        saveas(h, fullfile(datapath{i}, ['FatIntensityAll_' batchLabel{i} ext]));
    end
    
    % and the combined distributions
    clf
    for row = 1:nrows
        plotInfo = struct('plotTitle', plotTitles{row});
        subplotDist(nrows, ncols, row, plotInfo, statsTot, fieldnames{row}, batchLabel)
    end
    saveas(h, fullfile(combinedOutputPath, 'analysis_results', ['FatIntensityCombined' ext]));
    close
end

%%
%--------------------------------------------------------------
% reproduce visualization of full segmentation
%--------------------------------------------------------------
% 
% figure,
% 
% tic
% DsL = cellLayer.cellL(t, DsIdx);
% FatL = cellLayer.cellL(t, FatIdx);
% otherL = cellLayer.cellL(t, otherIdx);
% toc
% 
% % boundaries
% bondStates = cat(1,cellLayer.bonds{t}.state);
% bdryIdx = bondStates(:,4) == 1;
% tic
% bondL = cellLayer.bondL(t, bdryIdx);
% toc
% bdryMask = imdilate(bondL > 0, strel('disk',3));
% 
% % combine
% fullseg = cat(3,mat2gray(DsL>0 | bdryMask), FatL>0 | bdryMask, otherL>0 & ~bdryMask);
% 
% imshow(fullseg)

%%
%--------------------------------------------------------------
% interface counting
%--------------------------------------------------------------
t=1;
for i = goodseg

    cellLayer = CL{i};
    
    % store stats in output structure
    colorCount = sum(cat(1,cellLayer.cells{1}.state));

    stats{i}.NFatCells = colorCount(5);
    stats{i}.NDsCells = colorCount(6);

    % count fraction of Fat-Ds interfaces with boundaries
    bondStates = cat(1,cellLayer.bonds{t}.state);

    % FRACTIONAL VALUES, SOMETHING WENT WRONG?
    stats{i}.NFatDsInterfaces = round(sum(bondStates(:,1) == 1)/2);
    stats{i}.NFatFatInterfaces = round(sum(bondStates(:,2) == 1)/2);
    stats{i}.NDsDsInterfaces = round(sum(bondStates(:,3) == 1)/2);
    stats{i}.NFatDsBdries = round(sum(bondStates(:,4) == 1)/2);

    % fraction of FatDs interfaces that form a boundary
    bdryFraction = stats{i}.NFatDsBdries/stats{i}.NFatDsInterfaces;

    % fraction of Fat-Ds interfaces
    N = stats{i}.NDsDsInterfaces + stats{i}.NFatFatInterfaces + stats{i}.NFatDsInterfaces;

    [stats{i}.NDsDsInterfaces   stats{i}.NFatDsInterfaces;...
    stats{i}.NFatDsInterfaces  stats{i}.NFatFatInterfaces]./N;

    % expectated fractions for well mixed
    [colorCount(2)^2  2*colorCount(1)*colorCount(2);...
        2*colorCount(1)*colorCount(2) colorCount(1)^2]/sum(colorCount)^2;
end

%%
%--------------------------------------------------------------
% boundary fraction
%--------------------------------------------------------------

h = figure;
hold on
cmap = lines(nBatches);
batchAvg = zeros([nBatches 1]);
for i = 1:nBatches
    batchAvg(i) = 0;
    for j = 1:numel(batchIdx{i})
        NBd = stats{batchIdx{i}(j)}.NFatDsBdries;
        NIf = stats{batchIdx{i}(j)}.NFatDsInterfaces;
        
        if i == 6
            [j NBd/NIf]
        end
        scatter(batchTimes(i),NBd/NIf, 40, i, 'o', 'fill');
        batchAvg(i) = batchAvg(i) + NBd/NIf;
    end
    batchAvg(i) = batchAvg(i)/numel(batchIdx{i});
end
t = -2:24;
tau = 7;
minF = batchAvg(1);
maxF = batchAvg(nBatches-1);
lw = 2;
h = plot(t, minF + (maxF-minF)./(1+ exp(-(t-tau)/2)), 'k', 'Linewidth', lw);
hold off
%legend(h,'sigmoid function');
axis([t(1) t(end) 0 maxF + 0.1]);
ylabel('N-Fat-Ds / N-accumulation', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('time after Dox addition (hours)', 'FontSize', 14, 'FontWeight', 'bold');
title('Fraction of accumulating Fat-Ds interfaces', 'FontSize', 14, 'FontWeight', 'bold');
set(gcf, 'Color', 'w');
set(gca,'FontSize', 14)
set(gca,'FontWeight', 'bold')
set(gca, 'LineWidth', 2);

if saveResults
    %I = getframe(gcf);
    saveas(h, fullfile(combinedOutputPath, 'analysis_results', ['accumulationFraction' ext]));
%     for i = 1:nImages % for convenience make copy in each dir
%         saveas(h, fullfile(filepath{i}, 'analysis_results', 'accumulationFraction.png'));
%         %imwrite(I.cdata, fullfile(filepath{i}, 'analysis_results', 'DsIntensity.png'));
%     end
end

%%
%--------------------------------------------------------------
% total Ds and Fat level, density?
%--------------------------------------------------------------

h = figure,
hold on

cmap = lines(nBatches);
batchAvg = zeros([nBatches 1]);

for i = 1:nBatches
    batchAvg(i) = 0;
    for j = 1:numel(batchIdx{i})
        DsImean = mean(stats{batchIdx{i}(j)}.IDsSeg.table);
        scatter(batchTimes(i),DsImean, 40, i, 'o', 'fill');
        batchAvg(i) = batchAvg(i) + DsImean;
    end
    batchAvg(i) = batchAvg(i)/numel(batchIdx{i});
end

t = -2:24;
%tau = 6.5;
minDsI = batchAvg(1);
maxDsI = batchAvg(nBatches-2);
maxDsIf = batchAvg(nBatches);
lw = 2;
tau2 = 2;
%h = plot(t, minDsI + (maxDsI-minDsI)./(1+ exp(-(t-tau))),'k',...
%        t,  minDsI + (maxDsIf-minDsI)*(t-tau2)/(t(end)-tau2), '--k', 'Linewidth', lw);
hold off
%legend(h,{'sigmoid function','linear increase'}, 'Location', 'SouthEast');
axis([t(1) t(end) minDsI-50 batchAvg(nBatches-1)+100]);
ylabel('Ds level', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('time after Dox addition (hours)', 'FontSize', 14, 'FontWeight', 'bold');
title('Total Ds level in Ds Cells', 'FontSize', 14, 'FontWeight', 'bold');
set(gcf, 'Color', 'w');
set(gca,'FontSize', 14)
set(gca,'FontWeight', 'bold')
set(gca, 'LineWidth', 2);

if saveResults
    %I = getframe(gcf);
    saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['totalDsLevels' ext]));
%     for i = 1:nImages % for convenience make copy in each dir
%         saveas(h, fullfile(filepath{i}, 'analysis_results', 'accumulationFraction.png'));
%         %imwrite(I.cdata, fullfile(filepath{i}, 'analysis_results', 'DsIntensity.png'));
%     end
end

%%
%--------------------------------------------------------------
% bdries vs total Ds 
%--------------------------------------------------------------
figure
hold on

cmap = lines(nBatches);
batchAvg = zeros([nBatches 1]);

DsIall = [];
NBall = [];
NAcall = [];
statsPerTimeForDavid = {};

for i = 1:nBatches
    batchAvg(i) = 0;
    
    DsIbatch = [];
    NBbatch = [];
    NAcbatch = [];
    NDsCellsBatch = [];
    NFatCellsBatch = [];

    for j = 1:numel(batchIdx{i})
        
        NBd = stats{batchIdx{i}(j)}.NFatDsBdries;
        NIf = stats{batchIdx{i}(j)}.NFatDsInterfaces;
        DsImean = mean(stats{batchIdx{i}(j)}.IDsSeg.table);
        
        NFatCells = stats{batchIdx{i}(j)}.NFatCells;
        NDsCells = stats{batchIdx{i}(j)}.NFatCells;
        
        DsIall = [DsIall, DsImean];
        NBall = [NBall, NIf];
        NAcall = [NAcall, NBd];
        
        DsIbatch = [DsIbatch, DsImean];
        NBbatch = [NBbatch, NIf];
        NAcbatch = [NAcbatch, NBd];
        NFatCellsBatch = [NFatCellsBatch NFatCells];
        NDsCellsBatch = [NDsCellsBatch NDsCells];

        scatter(DsImean-minDsI,NBd/NIf, 40, i, 'o', 'fill');
    end
    
    statsPerTimeForDavid{i} = struct('label', batchLabel{i},...
                            'DsIntensities', DsIbatch,...
                            'NumBoundaries', NBbatch,...
                            'NumAccumulating', NAcbatch,...
                            'NFatCells',NFatCellsBatch,...
                            'NDsCells',NDsCellsBatch);
end
hold off

save(fullfile(combinedOutputPath,'bdryFractionDsIntensitiesForDavid'), 'statsPerTimeForDavid');

%axis([0 maxDsI-minDsI 0 maxF]);
ylabel('boundary fraction', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Ds level in Ds cells', 'FontSize', 14, 'FontWeight', 'bold');
title('Total Ds level vs. boundary fraction', 'FontSize', 14, 'FontWeight', 'bold');
set(gcf, 'Color', 'w');
set(gca,'FontSize', 14)
set(gca,'FontWeight', 'bold')
set(gca, 'LineWidth', 2);

if saveResults
    saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['DsVsBoundaries' ext]));
end

% %%
% load('/Users/idse/Dropbox/Sprinzak/shared/snapshots/snapshots for analysis/2.7.15/bdryFractionData.mat');
% scatter(S.DsIntensities, S.NumBoundaries./S.NumAccumulating,'.');
% axis([min(S.DsIntensities)  max(S.DsIntensities) 0 0.4]);

%%
%-----------------------------------------------------------------
% scatter plot and 2D histogram of Fat/Ds of cell in Fat/Ds interface and boundary
%-----------------------------------------------------------------

% for 2D histogram: 
% need intensities of cells on each side of bdry
% and also intensity of cells on each side of interface, but only if
% interface doesn't belong to cell that also has a boundary?

c1FatTot = [];
c2DsTot = [];
c1FatBTot = [];
c2DsBTot = [];

for i = 12%:nBatches

    c1FatbatchTot = [];
    c2DsbatchTot = [];
    c1FatBbatchTot = [];
    c2DsBbatchTot = [];
    
    for j = batchIdx{i}

        % background subtraction: minimal value in that image
        DsBG = 120; %min(stats{i}.IDs.table);
        FatBG = 110; %min(stats{i}.IFat.table);
        
        cellLayer = CL{j};

        % indices for bonds and interfaces
        bondStates = cat(1,cellLayer.bonds{1}.state);
        FatDsIfIdx = bondStates(:,1) == 1;
        FatDsBdryIdx = bondStates(:,4) == 1;

        ifCells = cat(1,cellLayer.bonds{1}(FatDsIfIdx).cellInd);
        bdryCells = cat(1,cellLayer.bonds{1}(FatDsBdryIdx).cellInd);
        
        if numel(bdryCells) > 1
            
            % every bond is there in two directions, so each pair is there twice
            % we throw the redundancy out by keeping the pairs where the first is
            % Fat
            c1state = cat(1,cellLayer.cells{1}(ifCells(:,1)).state);
            ifCells = ifCells(c1state(:,5) == 1,:);

            c1state = cat(1,cellLayer.cells{1}(bdryCells(:,1)).state);
            bdryCells = bdryCells(c1state(:,5) == 1,:);

            % all interfaces
            idx = ifCells;

            c1state = cat(1,cellLayer.cells{1}(idx(:,1)).state);
            c1Fat =  c1state(:,1);% - FatBG;

            c2state = cat(1,cellLayer.cells{1}(idx(:,2)).state);
            c2Ds =  c2state(:,2);% - DsBG;

            H = hist2d(real([c2Ds, c1Fat]), 50, 50, [Ibins{1}(1) Ibins{1}(end)], [Ibins{2}(1) Ibins{2}(end)]);
            H = H./sum(H(:));
            H = imrotate(H,90);
            loglims = [2 7.5];
            Hlog = hist2d(real([log(c2Ds), log(c1Fat)]), 50, 50, loglims, loglims);
            Hlog = Hlog./sum(Hlog(:));
            Hlog = imrotate(Hlog,90);

            % interfaces with accumulation
            idx = bdryCells;
        
            c1state = cat(1,cellLayer.cells{1}(idx(:,1)).state);
            c1FatB =  c1state(:,1);% - FatBG;
            c2state = cat(1,cellLayer.cells{1}(idx(:,2)).state);
            c2DsB =  c2state(:,2);% - DsBG;
            HB = hist2d(real([c2DsB, c1FatB]), 50, 50, [Ibins{1}(1) Ibins{1}(end)], [Ibins{2}(1) Ibins{2}(end)]);
            HB = HB./sum(HB(:));
            HB = imrotate(HB,90);
            
            HBlog = hist2d(real([log(c2DsB), log(c1FatB)]), 50, 50, loglims, loglims);
            HBlog = HBlog./sum(HBlog(:));
            HBlog = imrotate(HBlog,90);

            % collect batch total
            c1FatbatchTot = cat(1, c1FatbatchTot, c1Fat);
            c2DsbatchTot = cat(1, c2DsbatchTot, c2Ds);
            c1FatBbatchTot = cat(1, c1FatBbatchTot, c1FatB);
            c2DsBbatchTot = cat(1, c2DsBbatchTot, c2DsB);
        
%             % visualize
%             h = figure;
%             scatter(c1FatB,c2DsB,'oy', 'fill') 
%             hold on
%             scatter(c1Fat,c2Ds, 'ob') 
%             hold off
% 
%             title('Fat-Ds Pair Intensities');
%             xlabel('Fat');
%             ylabel('Ds');
%             axis([lims{1,2}, lims{1,1}]);
%             legend( ['Fat-Ds boundary, N = ' num2str(sum(FatDsBdryIdx)/2)],...
%                     ['Fat-Ds interface, N = ' num2str(sum(FatDsIfIdx)/2)]);

            if saveResults
                saveas(h, fullfile(filepath{j}, 'analysis_results', ['pairIntensities' ext]));
                clf
                imshow(4*cat(3,HB+H,HB,H)/max(H(:)));
                saveas(h, fullfile(filepath{j}, 'analysis_results', ['pairI2Dhist' ext]));

                imshow(4*cat(3,HBlog+Hlog,HBlog,Hlog)/max(Hlog(:)));
                saveas(h, fullfile(filepath{j}, 'analysis_results', ['pairI2DhistLog' ext]));
                %I = getframe(gcf);
                %imwrite(I.cdata, fullfile(filepath{i}, 'analysis_results', 'pairIntensities.png'));
            end
            close
            
        else
            warning([flabel{j} 'has one or no boundaries?']);
        end
    end

    %total total
    c1FatTot = cat(1, c1FatTot, c1FatbatchTot);
    c2DsTot = cat(1, c2DsTot, c2DsbatchTot);
	c1FatBTot = cat(1, c1FatBTot, c1FatBbatchTot);
    c2DsBTot = cat(1, c2DsBTot, c2DsBbatchTot);
    
	% batch hist
    HbatchTot = hist2d([c2DsbatchTot, c1FatbatchTot], 50, 50, [Ibins{1}(1) Ibins{1}(end)], [Ibins{2}(1) Ibins{2}(end)]);
    HbatchTot = HbatchTot./sum(HbatchTot(:));
    HbatchTot = imrotate(HbatchTot,90);
    
    HbatchTotlog = hist2d(real([log(c2DsTot), log(c1FatTot)]), 50, 50, loglims, loglims);
    HbatchTotlog = HbatchTotlog./sum(HbatchTotlog(:));
    HbatchTotlog = imrotate(HbatchTotlog,90);
    
    HBbatchTot = hist2d([c2DsBbatchTot, c1FatBbatchTot], 50, 50, [Ibins{1}(1) Ibins{1}(end)], [Ibins{2}(1) Ibins{2}(end)]);
    HBbatchTot = HBbatchTot./sum(HBbatchTot(:));
    HBbatchTot = imrotate(HBbatchTot,90);
    
    HBbatchTotlog = hist2d([log(c2DsBbatchTot), log(c1FatBbatchTot)], 50, 50, loglims, loglims);
    HBbatchTotlog = HBbatchTotlog./sum(HBbatchTotlog(:));
    HBbatchTotlog = imrotate(HBbatchTotlog,90);
    
    % visualize
    h = figure;
    scatter(c1FatBbatchTot,c2DsBbatchTot,'oy', 'fill') 
    hold on
    scatter(c1FatbatchTot,c2DsbatchTot, 'ob') 
    hold off
    title('Fat-Ds Pair Intensities');
    xlabel('Fat');
    ylabel('Ds');
    axis([lims{1,2}, lims{1,1}]);
    legend( 'Fat-Ds boundary','Fat-Ds interface');

    if saveResults
        
        saveas(h, fullfile(datapath{i},['pairIntensities' ext]));
        
        % yellow on top
        clf
        scatter(c1FatbatchTot,c2DsbatchTot, 'ob') 
        hold on
        scatter(c1FatBbatchTot,c2DsBbatchTot,'oy', 'fill') 
        hold off
        title('Fat-Ds Pair Intensities');
        xlabel('Fat');
        ylabel('Ds');
        axis([lims{1,2}, lims{1,1}]);
        legend( 'Fat-Ds boundary','Fat-Ds interface');
        saveas(h, fullfile(datapath{i},['pairIntensities2' ext]));
        
        % log
        clf
        scatter(log(c1FatBbatchTot),log(c2DsBbatchTot),'oy', 'fill') 
        hold on
        scatter(log(c1FatbatchTot),log(c2DsbatchTot), 'ob') 
        hold off
        axis([loglims loglims]);
        saveas(h, fullfile(datapath{i},['pairIntensitiesLog' ext]));
        
        % histogram
        clf
        HB = HBbatchTot;
        H = HbatchTot;
        imshow(4*cat(3,HB+H,HB,H)/max(H(:)));
        saveas(h, fullfile(datapath{i}, ['pairI2Dhist' ext]));
        imshow(mat2gray(H));
        saveas(h, fullfile(datapath{i}, ['pairI2DhistAll' ext]));
        imshow(mat2gray(HB));
        saveas(h, fullfile(datapath{i}, ['pairI2DhistBdry' ext]));

        % histogram log
        clf
        HB = HBbatchTotlog;
        H = HbatchTotlog;
        imshow(cat(3,HB+H,HB,H)/max(H(:)));
        saveas(h, fullfile(datapath{i}, ['pairI2DhistLog' ext]));
        imshow(mat2gray(H));

        saveas(h, fullfile(datapath{i}, ['pairI2DhistAllLog' ext]));
        imshow(mat2gray(HB));
        saveas(h, fullfile(datapath{i}, ['pairI2DhistBdryLog' ext]));
        %I = getframe(gcf);
        %imwrite(I.cdata, fullfile(filepath{i}, 'analysis_results', 'pairIntensities.png'));
    end
    close
end

% 
% %%
% %--------------------------------------------------------------
% % scatter plot Fat/Ds of boundary
% %--------------------------------------------------------------
% 
% for i = goodseg
% 
%     h = figure;
%     scatter(bdryI{i}.bdryFI, bdryI{i}.bdryDI)
%     
%     title('Fat-Ds Boundary Intensities');
%     xlabel('Fat');
%     ylabel('Ds');
%     axis([  lims{1,2}(1), 2*lims{1,2}(2),...
%             lims{1,1}(1), 2*lims{1,1}(2)]);
% 
%     if saveResults
%         saveas(h, fullfile(filepath{i}, 'analysis_results', ['bdryIntensities' ext]));
%         %I = getframe(gcf);
%         %imwrite(I.cdata, fullfile(filepath{i}, 'analysis_results', 'bdryIntensities.png'));
%     end
%     close
% end

% %% tracking down half integer bonds?
% 
% bondCell = cat(1,cellLayer.bonds{1}.cellInd);
% % bonds that represent cell interfaces should appear with both orientations
% % and thus be an even number:
% sum(bondCell(:,2)~=0)
% 
% % number of bonds that are Fat-Ds interfaces should therefore also be even
% bondStates = cat(1,cellLayer.bonds{t}.state);
% sum(bondStates(:,1)==1)
% 
% % but it is not...
% 
% %%
% %--------------------------------------------------------------
% % bond length longer for Fat-Ds interfaces?
% % also check cell geometry features
% %--------------------------------------------------------------
% 
% for i = goodseg %1:nImages
% 
%     cellLayer = CL{i};
%     
%     % bond geometry
%     %-----------------------------------
%     
%     cellLayer.calcBondGeometry(t);
%     bondLength = cellLayer.getBondGeometry(t, 'length');
%     bondStates = cat(1,cellLayer.bonds{t}.state);
%     
%     % the others are on im edge
%     %realBonds = cat(1,cellLayer.bonds{1}.label)~=0; 
%     cellInd = cat(1,cellLayer.bonds{1}.cellInd);
%     interfaceBonds = cellInd(:,2) ~= 0;
%     L = bondLength(interfaceBonds);
%     L = L(~isnan(L)); % FIGURE OUT WHY THIS IS NECESSARY
% 
%     FDL = bondLength(bondStates(:,1) == 1,:);
%     FDL = FDL(~isnan(FDL));
% 
%     FFL = bondLength(bondStates(:,2) == 1,:);
%     FFL = FFL(~isnan(FFL));
%     
%     DDL = bondLength(bondStates(:,3) == 1,:);
%     DDL = DDL(~isnan(DDL));
%     
%     BL = bondLength(bondStates(:,4) == 1,:);
%     BL = BL(~isnan(BL));
% 
%     stats{i}.ifLength = makeDist(L, 20);
%     stats{i}.FDifLength = makeDist(FDL, stats{i}.ifLength.bins);
%     stats{i}.FFifLength = makeDist(FFL, stats{i}.ifLength.bins);
%     stats{i}.DDifLength = makeDist(DDL, stats{i}.ifLength.bins);
%     stats{i}.FDbdryLength = makeDist(BL, stats{i}.ifLength.bins);
% 
%     % visualize
%     ncols = 2;
%     nrows = 3;
% 
%     scrsz = get(0,'ScreenSize');
%     h = figure;%('Position',[1 scrsz(4)/2 scrsz(3)*0.7 scrsz(4)])
% 
%     plotTitles = {'interface length', 'Fat-Ds if length', 'Fat-Fat if length','Ds-Ds if length', 'Fat-Ds bdry length'};
%     fieldnames = {'ifLength','FDifLength','FFifLength','DDifLength','FDbdryLength'};
%     flabel = {''};
% 
%     row = 1;
%     plotInfo = struct('plotTitle', plotTitles{row});
%     subplotDistCompare(nrows, ncols, row, plotInfo, stats, fieldnames, i)
%     
%     % cell geometry
%     %-----------------------------------
%     
%     cellLayer.calcCellGeometry(t);
%     cA = cellLayer.getCellGeometry(t, 'area');
%     majAx = cellLayer.getCellGeometry(t, 'majorAxisLength');
%     minAx = cellLayer.getCellGeometry(t, 'minorAxisLength');
%     anisotropy = (majAx - minAx)./(majAx+minAx);
%     
%     % separate by type
%     cellStates = cat(1,cellLayer.cells{t}.state);
%     FatIdx = cellStates(:,5)==1;
%     DsIdx = cellStates(:,6)==1;
%     otherIdx = cellStates(:,5)~=1 & cellStates(:,6)~=1;
%     
%     stats{i}.area = makeDist(cA, 20);
%     stats{i}.Farea = makeDist(cA(FatIdx), stats{i}.area.bins);
%     stats{i}.Darea = makeDist(cA(DsIdx), stats{i}.area.bins);
% 
%     stats{i}.anisotropy = makeDist(anisotropy, 20);
%     stats{i}.Fanisotropy = makeDist(anisotropy(FatIdx), stats{i}.anisotropy.bins);
%     stats{i}.Danisotropy = makeDist(anisotropy(DsIdx), stats{i}.anisotropy.bins);
%     
%     % visualize
%     plotTitles = {'cell area', 'Fat area', 'Ds area'};
%     fieldnames = {'area', 'Farea', 'Darea'};
%     flabel = {''};
% 
%     row = 2;
%     plotInfo = struct('plotTitle', plotTitles{row});
%     subplotDistCompare(nrows, ncols, row, plotInfo, stats, fieldnames, i)
%     
%     plotTitles = {'cell anisotropy', 'Fat anisotropy', 'Ds anisotropy'};
%     fieldnames = {'anisotropy', 'Fanisotropy', 'Danisotropy'};
%     flabel = {''};
% 
%     row = 3;
%     plotInfo = struct('plotTitle', plotTitles{row});
%     subplotDistCompare(nrows, ncols, row, plotInfo, stats, fieldnames, i)
%     
%     % save the whole thing
%     set(gcf, 'Color', 'w');
%     if saveResults
%         I = getframe(gcf);
%         saveas(h, fullfile(filepath{i}, 'analysis_results', 'ifLength.png'));
%         %imwrite(I.cdata, fullfile(filepath{i}, 'analysis_results', 'ifLength.png'));
%     end
%     close
% end
% %mmdKS(stats{i}.ifLength.table, stats{i}.FDifLength.table)
% 
% %%
% %--------------------------------------------------------------
% % get the radial distribution functions
% %--------------------------------------------------------------
% 
% for i = 1:nImages
%     
%     cellLayer = CL{i};
% 
%     if isempty(cellLayer.cells{t}(1).geometry)
%         cellLayer.calcCellGeometry(t);
%     end
% 
%     centroids = cellLayer.getCellGeometry(t, 'centroid')';
% 
%     % cell type indices
%     cellStates = cat(1,cellLayer.cells{1}.state);
%     FatIdx = cellStates(:,5)==1;
%     DsIdx = cellStates(:,6)==1;
% 
%     % overall density
%     yrange = [0 size(cellLayer.L,1)];
%     xrange = [0 size(cellLayer.L,2)];
% 
%     % radial binning
%     dr = 10;
%     rmax = 140;
%     rbins = 0:dr:rmax;
% 
%     % cell positions
%     X = centroids;
%     XF = centroids(FatIdx,:);
%     XD = centroids(DsIdx,:);
% 
%     % the distribution functions
%     g = radialDistributionFunction(X,X,xrange,yrange,rbins);
%     gFF = radialDistributionFunction(XF,XF,xrange,yrange,rbins);
%     gDD = radialDistributionFunction(XD,XD,xrange,yrange,rbins);
%     gFD = radialDistributionFunction(XF,XD,xrange,yrange,rbins);
%     gDF = radialDistributionFunction(XD,XF,xrange,yrange,rbins);
% 
%     % visualize
%     h = figure;
%     plot(rbins(1:end-1)+dr/2, g, 'b', 'LineWidth', 2);
%     hold on
%     plot(rbins(1:end-1)+dr/2,gFF, 'g', 'LineWidth', 2)
%     plot(rbins(1:end-1)+dr/2,gDD, 'r', 'LineWidth', 2)
%     plot(rbins(1:end-1)+dr/2,(gFD+gDF)/2, 'm', 'LineWidth', 2)
%     hold off
%     legend('g','gFF','gDD','gFD');
%     
%     % save the whole thing
%     set(gcf, 'Color', 'w');
%     if saveResults
%         saveas(h, fullfile(filepath{i}, 'analysis_results', 'radialDistFn.png'));
% %         I = getframe(gcf);
% %         imwrite(I.cdata, fullfile(filepath{i}, 'analysis_results', 'radialDistFn.png'));
%     end
%     close
% end

%%
%-----------------------------------------------------------------
% scatter plot of Ds of cell vs boundary intensity
%-----------------------------------------------------------------

c1FatTot = [];
c2DsTot = [];
c1FatBTot = [];
c2DsBTot = [];
bdryDsTot = [];

for i = 5%:nBatches

    c1FatbatchTot = [];
    c2DsbatchTot = [];
    c1FatBbatchTot = [];
    c2DsBbatchTot = [];
    bdryDsbatchTot = [];
    
    for j = batchIdx{i}

        % background subtraction: minimal value in that image
        DsBG = 120; %min(stats{i}.IDs.table);
        FatBG = 110; %min(stats{i}.IFat.table);
        
        cellLayer = CL{j};

        % indices for bonds and interfaces
        bondStates = cat(1,cellLayer.bonds{1}.state);
        FatDsIfIdx = bondStates(:,1) == 1;
        FatDsBdryIdx = bondStates(:,4) == 1;

        ifCells = cat(1,cellLayer.bonds{1}(FatDsIfIdx).cellInd);
        bdryCells = cat(1,cellLayer.bonds{1}(FatDsBdryIdx).cellInd);
        bdryState = cat(1,cellLayer.bonds{1}(FatDsBdryIdx).state);
        bdryDsI = bdryState(:,6);
        
        if numel(bdryCells) > 1
            
            % every bond is there in two directions, so each pair is there twice
            % we throw the redundancy out by keeping the pairs where the first is
            % Fat
            c1state = cat(1,cellLayer.cells{1}(ifCells(:,1)).state);
            ifCells = ifCells(c1state(:,5) == 1,:);

            c1state = cat(1,cellLayer.cells{1}(bdryCells(:,1)).state);
            bdryCells = bdryCells(c1state(:,5) == 1,:);

            bdryDsI = bdryDsI(c1state(:,5) == 1);
            
            % interfaces with accumulation
            idx = bdryCells;

            c1state = cat(1,cellLayer.cells{1}(idx(:,1)).state);
            c1FatB =  c1state(:,1) - FatBG;
            c2state = cat(1,cellLayer.cells{1}(idx(:,2)).state);
            c2DsB =  c2state(:,2) - DsBG;

            % collect batch total
            c1FatBbatchTot = cat(1, c1FatBbatchTot, c1FatB);
            c2DsBbatchTot = cat(1, c2DsBbatchTot, c2DsB);
            bdryDsbatchTot = cat(1, bdryDsbatchTot, bdryDsI); 
            
            % visualize
            h = figure;
            scatter(c1FatB,c2DsB,'o', 'fill') 
%             hold on
%             scatter(c1Fat,c2Ds, 'ob') 
%             hold off

            title('Fat-Ds Pair Intensities');
            xlabel('Fat');
            ylabel('Ds');
            axis([lims{1,2}, lims{1,1}]);
            legend( ['Fat-Ds boundary, N = ' num2str(sum(FatDsBdryIdx)/2)],...
                    ['Fat-Ds interface, N = ' num2str(sum(FatDsIfIdx)/2)]);
            close
        else
            warning([flabel{j} 'has one or no boundaries?']);
        end
    end

    %total total
	c1FatBTot = cat(1, c1FatBTot, c1FatBbatchTot);
    c2DsBTot = cat(1, c2DsBTot, c2DsBbatchTot);
    bdryDsTot = cat(1, bdryDsTot, bdryDsbatchTot);
    
    % visualize
    h = figure;
    scatter(c1FatBbatchTot,bdryDsTot,'.') 
    title('cell Ds vs boundary Ds Intensities');
    xlabel('cell Ds');
    ylabel('boundary Ds');
    corr2(c1FatBbatchTot,bdryDsTot)
    axis([lims{1,2}, lims{1,1}]);
    %legend( 'Fat-Ds boundary','Fat-Ds interface');
end

%%
%------------------------------------------------
% save stats human readable and matlab format
%------------------------------------------------

for i = 1:nImages

    txtfile = fullfile(filepath{i}, 'analysis_results', 'stats_new.txt');
    if exist(txtfile)
        delete(fullfile(filepath{i}, 'analysis_results', 'stats_new.txt'));
    end
    diary(txtfile)
%     i
%     stats{i}
    diary off
    S = stats{i};
    save(fullfile(filepath{i}, 'analysis_results', 'stats_new'), 'S');
end

%%
%-------------------------------------------------------
% save stats human readable and matlab format for elife
%-------------------------------------------------------

for i = 1:nImages

    fparts = strsplit(filepath{i},'/');
    fname = fullfile(combinedOutputPath, ['stats_' fparts{end-1} '_' fparts{end}]);
    txtfile = [fname '.txt'];

    if exist(txtfile)
        delete(txtfile);
    end
    diary(txtfile)
    S = stats{i};
    fn_structdisp(S)
    diary off

%     S = stats{i};
%     save(fname, 'S');
end

%%
%-------------------------------------------------------
% save stats in single .txt for elife
%-------------------------------------------------------

fname = fullfile(combinedOutputPath, 'allstats');
txtfile = [fname '.txt'];

if exist(txtfile)
    delete(txtfile);
end
diary(txtfile)

for i = 1:nImages

    fparts = strsplit(filepath{i},'/');
    disp('----------------');
    disp(fparts{end})
    disp('----------------');
    S = stats{i};
    fn_structdisp(S)
end

diary off