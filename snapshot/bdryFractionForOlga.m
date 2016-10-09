clear all;

load(fullfile('/Users/idse/','bdryFrac.mat'));
saveResults = false;

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
%h = plot(t, minF + (maxF-minF)./(1+ exp(-(t-tau)/2)), 'k', 'Linewidth', lw);
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
% total Ds level
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
