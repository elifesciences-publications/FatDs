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

for i = 1:nImages
    fname = fullfile(filepath{i},'matlab_seg',[flabel{i} '_seg617']); %_seg_new
    if exist([fname '.mat'],'file')
        S = load(fname);
        CL{i} = S.cellLayer;
    else
        warning([fname ' doesnt exist']);
    end
end

% save results
%--------------------------------------------------
saveResults = false;

% PAY ATTENTION: there was _new at the end here at some point
ext = '526.png';

stats = {};

% fixed samples, not time
ti = 1;

Nbins = 50;

%%
%--------------------------------------------------------------
% collect statistics
%--------------------------------------------------------------

stats = {};
statsTot = {};

for i = 1:nBatches
    
    IDsTot = [];
    IFatTot = [];
    
    IDsBdryTot = [];
    IFatBdryTot = [];
    IDsNoBdryTot = [];
    IFatNoBdryTot = [];

    IDsSegTot = [];
    IFatSegTot = [];
    
    ifIFatDsTot = [];
    bdryIFatDsTot = [];
    
    for j = batchIdx{i}

        cellLayer = CL{j};

        % cell and bond states
        cellStates = cat(1,cellLayer.cells{ti}.state);
        bondStates = cat(1,cellLayer.bonds{ti}.state);

        % indices for subsets of cells and boundaries
        %------------------------------------------------------
        
        FatIdx = cellStates(:,5)==1;
        DsIdx = cellStates(:,6)==1;

        FatDsIfIdx = bondStates(:,1) == 1;
        FatFatIfIdx = bondStates(:,2) == 1;
        DsDsIfIdx = bondStates(:,3) == 1;
        FatDsBdryIdx = bondStates(:,4) == 1;
        
        % cells with Fat-Ds interfaces or accumulating bdries
        ifCells = cat(1,cellLayer.bonds{ti}(FatDsIfIdx).cellInd);
        bdryCells = cat(1,cellLayer.bonds{ti}(FatDsBdryIdx).cellInd);
        
        % every bond is there in two directions, so each pair is there twice
        % we throw the redundancy out by keeping the pairs where the first
        % is Fat
        c1state = cat(1,cellLayer.cells{1}(ifCells(:,1)).state);
        ifCells = ifCells(c1state(:,5) == 1,:);
        c1state = cat(1,cellLayer.cells{1}(bdryCells(:,1)).state);
        bdryCells = bdryCells(c1state(:,5) == 1,:);

        rmidx = false([size(ifCells,1) 1]);
        for ifi = 1:size(ifCells,1)
            if any(bdryCells(:,1) == ifCells(ifi,1))
                rmidx(ifi) = true;
            end
        end
        ifonlyCells = ifCells;
        ifonlyCells(rmidx,:)=[];
        %[numel(ifCells) numel(ifonlyCells) sum(FatDsBdryIdx)/sum(FatDsIfIdx)]
        ifCells = ifonlyCells;
        
        % to avoid confusion: bdryCells are cell pairs belonging to
        % specific boundaries, cellWBdry, are cells that have a boundary
        % so the first may contain the same cell multiple times
        
        % index of cells that have a boundary
        % index of cells that have no boundary but do have a Fat-Ds interface
        cellsWBdry = false(size(FatIdx));
        cellsWFDif = false(size(FatIdx)); 
        
        for ci = 1:numel(cellsWBdry)

            % bond states in cell
            cbi = cellLayer.cells{ti}(ci).bondInd;
            bsOfCell = cat(1, cellLayer.bonds{ti}(cbi).state);

            % 4th column of bond.state is boolean bdry
            if ~isempty(bsOfCell) && any(bsOfCell(:,4)==1)
                cellsWBdry(ci)=true;
            end 
            if ~isempty(bsOfCell) && any(bsOfCell(:,1)==1)
                cellsWFDif(ci)=true;
            end 
        end
        %[sum(cellsWBdry) sum(cellsWFDif) sum(cellsWBdry)/sum(cellsWFDif)]

        % total numbers of cells and interfaces
        %------------------------------------------------------
        
        % count fraction of Fat-Ds interfaces with boundaries
        stats{j}.NFatCells = sum(FatIdx);
        stats{j}.NDsCells = sum(DsIdx);
        stats{j}.NFatDsInterfaces = round(sum(FatDsIfIdx)/2);
        stats{j}.NFatFatInterfaces = round(sum(FatFatIfIdx)/2);
        stats{j}.NDsDsInterfaces = round(sum(DsDsIfIdx)/2);
        stats{j}.NFatDsBdries = round(sum(FatDsBdryIdx)/2);
       
        % tables of cell intensities
        %------------------------------------------------------
        
        IFat = cellStates(:,1);
        IDs = cellStates(:,2);

        stats{j}.IDsSeg = IDs(DsIdx);
        stats{j}.IFatSeg = IFat(FatIdx);

        % BDRY == SEG & BDRY
        IDsBdryTot = cat(1, IDsBdryTot, IDs(cellsWBdry & DsIdx));
        IFatBdryTot = cat(1, IFatBdryTot, IFat(cellsWBdry & FatIdx)); 

        % NO BDRY == SEG & NO BDRY
        IDsNoBdryTot = cat(1, IDsNoBdryTot, IDs(~cellsWBdry & DsIdx));
        IFatNoBdryTot = cat(1, IFatNoBdryTot, IFat(~cellsWBdry & FatIdx)); 
        
        % ALL: intensity distribution of all cells
        IDsTot = cat(1, IDsTot, IDs);
        IFatTot = cat(1, IFatTot, IFat);

        % SEG: intensity distribution of segmented cells
        IDsSegTot = cat(1, IDsSegTot, IDs(DsIdx));
        IFatSegTot = cat(1, IFatSegTot, IFat(FatIdx));
        
        % intensities of cells belonging to Fat-Ds interfaces with and
        % without accumulation
        if numel(bdryCells) > 1

            % ifCells was already set up to have cell 1 be the Fat cell
            c1state = cat(1,cellLayer.cells{ti}(ifCells(:,1)).state);
            c2state = cat(1,cellLayer.cells{ti}(ifCells(:,2)).state);
            ifIFatDs =  [c1state(:,1) c2state(:,2)];
            ifIFatDsTot = cat(1, ifIFatDsTot, ifIFatDs);
            
            c1state = cat(1,cellLayer.cells{ti}(bdryCells(:,1)).state);
            c2state = cat(1,cellLayer.cells{ti}(bdryCells(:,2)).state);
            bdryIFatDs =  [c1state(:,1) c2state(:,2)];
            bdryIFatDsTot = cat(1, bdryIFatDsTot, bdryIFatDs);
        end
    end

    % set up bins for distribution
    %-----------------------------

    lambdaI = 10; 
    tol = 0.005;
    IDsmin = min(IDsSegTot) - lambdaI;       
    IDsmax = (max(IDsSegTot)-min(IDsSegTot))*stretchlim(mat2gray(IDsSegTot), tol);
    IDsmax = IDsmax(2);

    tol = 0.005;
    IFatmin = min(IFatSegTot) - lambdaI;       
    IFatmax = (max(IFatSegTot)- min(IFatSegTot))*stretchlim(mat2gray(IFatSegTot), tol);
    IFatmax = IFatmax(2);

    DsBins = linspace(lambdaI, IDsmax, Nbins);
    FatBins = linspace(lambdaI, IFatmax, Nbins);
    
    logbinsDs = linspace(log(lambdaI), log(IDsmax), Nbins);
    logbinsFat = linspace(log(lambdaI), log(IFatmax), Nbins);
    
    % check:
    % sum(IDsSegTot - IDsmin > IDsmax)/numel(IDsSegTot)

    % make combined distribution
    %-----------------------------
    
    batchStat = struct();
    
    % ALL
    batchStat.IDs = makeDist(IDsTot, DsBins);
    batchStat.IFat = makeDist(IFatTot, FatBins);
     
    % SEG
    batchStat.IDsSeg = makeDist(IDsSegTot, DsBins);
    batchStat.IFatSeg = makeDist(IFatSegTot, FatBins);

    % SEG NO BDRY
    batchStat.IDsNoBdry = makeDist(IDsNoBdryTot - IDsmin, DsBins);
    batchStat.IFatNoBdry = makeDist(IFatNoBdryTot - IFatmin, FatBins);

    % SEG NO BDRY LOG
    batchStat.IDsNoBdryLogBin = makeDist(log(IDsNoBdryTot - IDsmin), logbinsDs);
    batchStat.IFatNoBdryLogBin = makeDist(log(IFatNoBdryTot - IFatmin), logbinsFat);
          
    % SEG BDRY
    batchStat.IDsBdry = makeDist(IDsBdryTot - IDsmin, DsBins);
    batchStat.IFatBdry = makeDist(IFatBdryTot - IFatmin, FatBins);

    % SEG BDRY LOG
    batchStat.IDsBdryLogBin = makeDist(log(IDsBdryTot  - IDsmin), logbinsDs);
    batchStat.IFatBdryLogBin = makeDist(log(IFatBdryTot - IFatmin), logbinsFat);

    % 2D distributions for interfaces
    I = [ifIFatDsTot(:,1) - IFatmin, ifIFatDsTot(:,2) - IDsmin];
    batchStat.IifLogBinHist = hist2d(log(I)', Nbins, Nbins, logbinsFat([1 end]), logbinsDs([1 end]));
    
    I = [bdryIFatDsTot(:,1) - IFatmin, bdryIFatDsTot(:,2) - IDsmin];
    batchStat.IbdryLogBinHist = hist2d(log(I), Nbins, Nbins, logbinsFat([1 end]), logbinsDs([1 end]));
    
    statsTot{i} = batchStat;
    
    % SAVE intensity tables to csv
    csvwrite(fullfile(combinedOutputPath, 'analysis_results', 'raw', ['IDs_' batchLabel{i} '.csv']), IDsTot);
    csvwrite(fullfile(combinedOutputPath, 'analysis_results', 'raw', ['IFat_' batchLabel{i} '.csv']), IFatTot);
    csvwrite(fullfile(combinedOutputPath, 'analysis_results', 'raw', ['IFatDsAccumulating_' batchLabel{i} '.csv']), bdryIFatDsTot);
    csvwrite(fullfile(combinedOutputPath, 'analysis_results', 'raw', ['IFatDsNonAccumulating_' batchLabel{i} '.csv']), ifIFatDsTot);
end

%% check that distributions are made correctly and compare 2D dist

bi = 12;

figure,
% equivalence of binning logI linearly or I on log scale
dist = hist(log(IDsNoBdryTot - IDsmin), logbinsDs);
dist = dist./sum(dist(:));
semilogx(exp(logbinsDs), dist,'-x','Color','r')
xlim([10 DsBins(end)]);

figure,
dist = hist(IDsNoBdryTot - IDsmin, exp(logbinsDs));
dist = dist./sum(dist(:));
semilogx(exp(logbinsDs), dist,'-x')
xlim([10 DsBins(end)]);

% Log Bin
figure,
bins = exp(logbinsDs(1:end-1));
semilogx(bins, statsTot{bi}.IDsNoBdryLogBin.dist(1:end-1), '-x', 'Color', 'r', 'LineWidth', 2);
hold on
semilogx(bins, statsTot{bi}.IDsBdryLogBin.dist(1:end-1), '-x', 'Color', 'b', 'LineWidth', 2);
xlim([10 bins(end)]);
title('logbin');
hold off

% Linear
figure,
plot(DsBins(1:end-1), statsTot{bi}.IDsNoBdry.dist(1:end-1), '-x', 'Color', 'r', 'LineWidth', 2);
hold on
plot(DsBins(1:end-1), statsTot{bi}.IDsBdry.dist(1:end-1), '-x', 'Color', 'b', 'LineWidth', 2);
xlim([10 DsBins(end)]);
title('linear');
hold off

% Linear on Log Scale
figure,
semilogx(DsBins(1:end-1), statsTot{bi}.IDsNoBdry.dist(1:end-1), '-x', 'Color', 'r', 'LineWidth', 2);
hold on
semilogx(DsBins(1:end-1), statsTot{bi}.IDsBdry.dist(1:end-1), '-x', 'Color', 'b', 'LineWidth', 2);
xlim([10 DsBins(end)]);
title('logplot, regular bin');
hold off


%% look at Fat avg over time and Fat median

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

%% only 20 h combined distribution Ds for figure

figure,

bins = statsTot{bi}.IDsNoBdry.bins;
nobdrydist = statsTot{bi}.IDsNoBdry.dist;
bdrydist = statsTot{bi}.IDsBdry.dist;

[x,y] = histForBarlikePlot(bins, nobdrydist');
plot(x, y, 'Color', 'r', 'LineWidth', 2);
hold on
[x,y] = histForBarlikePlot(bins, bdrydist');
plot(x, y, 'Color', 'k', 'LineWidth', 2);
hold off
xlim([10 bins(end)])
legend({'Ds no boundary','Ds boundary'})
if saveResults 
    saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['DsIntensity20hAvgCombined' ext]));
end

%% only 20 h combined distribution Ds for figure

figure, 

bins = statsTot{bi}.IDsNoBdry.bins;
nobdrydist = statsTot{bi}.IDsNoBdry.dist;
bdrydist = statsTot{bi}.IDsBdry.dist;

[x,y] = histForBarlikePlot(bins, nobdrydist');
semilogx(x, y, 'Color', 'r', 'LineWidth', 2);
hold on
[x,y] = histForBarlikePlot(bins, bdrydist');
semilogx(x, y, 'Color', 'k', 'LineWidth', 2);
hold off
xlim([10 bins(end)])
legend({'Ds no boundary','Ds boundary'})
if saveResults 
    saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['DsIntensity20hAvgCombinedLog' ext]));
end

%% only 20 h combined distribution Ds on log scale for figure

figure,

for bi = 1:nBatches
    
    clf

    bins = exp(logbinsDs(1:end-1));
    nobdrydist = statsTot{bi}.IDsNoBdryLogBin.dist(1:end-1);
    bdrydist = statsTot{bi}.IDsBdryLogBin.dist(1:end-1);

    [x,y] = histForBarlikePlot(bins, nobdrydist');
    semilogx(x, y, '-', 'Color', 'r', 'LineWidth', 2);
    hold on
    [x,y] = histForBarlikePlot(bins, bdrydist');
    semilogx(x, y, '-', 'Color', 'k', 'LineWidth', 2);
    hold off
    legend({'Ds no boundary','Ds boundary'})
    axis([10 bins(end) 0 0.075])
    
    if saveResults
        saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['DsIntensity' batchLabel{bi} 'LogBinned.png']));
        saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['DsIntensity' batchLabel{bi} 'LogBinned.fig']));
    end
end

%% only 20 h combined distribution Fat

figure

bins = statsTot{bi}.IFatNoBdry.bins(1:end-1);
nobdrydist = statsTot{12}.IFatNoBdry.dist(1:end-1);
bdrydist = statsTot{12}.IFatBdry.dist(1:end-1);

[x,y] = histForBarlikePlot(bins, nobdrydist');
plot(x, y, '-', 'Color', 'r', 'LineWidth', 2);
hold on
[x,y] = histForBarlikePlot(bins, bdrydist');
plot(x, y, '-', 'Color', 'k', 'LineWidth', 2);
hold off
xlim([10 bins(end)]);
legend({'Fat no boundary','Fat boundary'})
if saveResults
    saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['FatIntensity20hAvgCombined' ext]));
end

%% only 20 h combined distribution Fat, log

figure 

bins = statsTot{bi}.IFatNoBdry.bins(1:end-1);
nobdrydist = statsTot{12}.IFatNoBdry.dist(1:end-1);
bdrydist = statsTot{12}.IFatBdry.dist(1:end-1);

[x,y] = histForBarlikePlot(bins, nobdrydist');
semilogx(x, y, '-', 'Color', 'r', 'LineWidth', 2);
hold on
[x,y] = histForBarlikePlot(bins, bdrydist');
semilogx(x, y, '-', 'Color', 'k', 'LineWidth', 2);
hold off
xlim([10 bins(end)]);
legend({'Fat no boundary','Fat boundary'})
if saveResults
    saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['FatIntensity20hAvgCombinedLog' ext]));
end

%% only 20 h combined distribution Fat, log binned

figure

for bi = 1:nBatches
    
    clf

    bins = exp(logbinsFat(1:end-1));
    nobdrydist = statsTot{bi}.IFatNoBdryLogBin.dist(1:end-1);
    bdrydist = statsTot{bi}.IFatBdryLogBin.dist(1:end-1);

    [x,y] = histForBarlikePlot(bins, nobdrydist');
    semilogx(x, y, '-', 'Color', 'r', 'LineWidth', 2);
    hold on
    [x,y] = histForBarlikePlot(bins, bdrydist');
    semilogx(x, y, '-', 'Color', 'k', 'LineWidth', 2);
    hold off
    axis([10 bins(end) 0 0.075])

    legend({'Fat no boundary','Fat boundary'})

    if saveResults
        saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['FatIntensity' batchLabel{bi} 'LogBinned.fig']));
        saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['FatIntensity' batchLabel{bi} 'LogBinned.png']));
    end
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
% bondStates = cat(1,cellLayer.bonds{ti}.state);
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
        
        if i == 12
            [j NBd/NIf];
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
    saveas(h, fullfile(combinedOutputPath, 'analysis_results', 'accumulationFraction.png'));
    saveas(h, fullfile(combinedOutputPath, 'analysis_results', 'accumulationFraction.fig'));
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
        DsImean = mean(stats{batchIdx{i}(j)}.IDsSeg);
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
    saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', 'totalDsLevels.png'));
    saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', 'totalDsLevels.fig'));
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
        DsImean = mean(stats{batchIdx{i}(j)}.IDsSeg);
        
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

save(fullfile(combinedOutputPath,'analysis_results','bdryFractionDsIntensitiesForDavid'), 'statsPerTimeForDavid');

%axis([0 maxDsI-minDsI 0 maxF]);
ylabel('boundary fraction', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Ds level in Ds cells', 'FontSize', 14, 'FontWeight', 'bold');
title('Total Ds level vs. boundary fraction', 'FontSize', 14, 'FontWeight', 'bold');
set(gcf, 'Color', 'w');
set(gca,'FontSize', 14)
set(gca,'FontWeight', 'bold')
set(gca, 'LineWidth', 2);

if saveResults
    saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', 'DsVsBoundaries.fig'));
    saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', 'DsVsBoundaries.png'));
end

% %%
% load('/Users/idse/Dropbox/Sprinzak/shared/snapshots/snapshots for analysis/2.7.15/bdryFractionData.mat');
% scatter(S.DsIntensities, S.NumBoundaries./S.NumAccumulating,'.');
% axis([min(S.DsIntensities)  max(S.DsIntensities) 0 0.4]);

%%
%-----------------------------------------------------------------
% 2D histogram of Fat/Ds of cell in Fat/Ds interface and boundary
%-----------------------------------------------------------------

for bi = 1:nBatches
        
        H = statsTot{bi}.IifLogBinHist;
        HB = statsTot{bi}.IbdryLogBinHist;
        
        T = (sum(HB(:))+sum(H(:)));
        H = H./T;
        HB = HB./T;

        % histogram log
        clf
%         tol = 0.02;
%         HBlim = stretchlim(mat2gray(HB),tol);
%         Hlim = stretchlim(mat2gray(H),tol);
%         HB = imadjust(mat2gray(HB, [10^(-4) 2*10^(-3)]))';%, HBlim)';
%         H = imadjust(mat2gray(H, [10^(-4) 2*10^(-3)]))';%, Hlim)';
        HB = mat2gray(HB, [10^(-4) 0.5*10^(-3)])';
        H = mat2gray(H, [10^(-4) 10^(-3)])';
   
        %iptsetpref('ImshowAxesVisible','off');
        imshow(cat(3,HB+H,HB,H));%,'XData',exp(logbinsDs))
        %axis square
        set(gca,'Ydir','normal');
        %set(gca,'Xscale','log','Xdir','normal');
        ylabel('Fat')
        xlabel('Ds')
        
        saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['2DhistLogV2' batchLabel{bi} '.fig']));
        saveas(gcf, fullfile(combinedOutputPath, 'analysis_results', ['2DhistLogV2' batchLabel{bi} '.png']));
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