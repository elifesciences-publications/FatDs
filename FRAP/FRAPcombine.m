%% plot of values from different movies

FRAPdataLocations
D = zeros([1 numel(movies)]);
sigD = zeros([1 numel(movies)]);
gamma = zeros([1 numel(movies)]);
sigGamma = zeros([1 numel(movies)]);
nobueno = D;

diary(fullfile(datapath,'FRAPfitValues.txt'));
for movieIdx = 1:numel(movies)
    
    disp(movies{movieIdx});
    if ~exist(fullfile(datapath, movies{movieIdx}),'file')
        disp('does not exist, moving on');
        nobueno(movieIdx) = true;
    else
        filepath = fileparts(fullfile(datapath, movies{movieIdx})); 
        load(fullfile(filepath,'result.mat'));
        
        D(movieIdx) = allParam.D(1);
        sigD(movieIdx) = allParam.D(2);
        disp(['D = ' num2str(D(movieIdx),'%.2e') '(' num2str(sigD(movieIdx),'%.2e') ')']);
        
        gamma(movieIdx)  = allParam.gamma(1);
        sigGamma(movieIdx)  = allParam.gamma(2);
        disp(['gamma = ' num2str(gamma(movieIdx),'%.2e') '(' num2str(sigGamma(movieIdx),'%.2e') ')']);
        %Draw(movieIdx) = allParam.Draw(1);
        %sigDraw(movieIdx) = allParam.Draw(2);
    end
    disp(' ');
end
diary off;

sigD = sigD(~nobueno);
D = D(~nobueno);
gamma = gamma(~nobueno);
sigGamma = sigGamma(~nobueno);

N = sum(~nobueno);

%% combine diffusion constants

w = 1./(sigD.^2);
w = w./sum(w);
Dtot = sum(w.*D);
sigTot = 1/sqrt(sum(1./(sigD.^2)));

h = errorbar(1:N, D, sigD,'.');
hold on
errorbar(N+1, Dtot, sigTot,'.r');
hold off
title(['D from different movies, Dtot = ' num2str(Dtot,'%.2e') '(' num2str(sigTot,'%.2e') ')']);
saveas(gcf, fullfile(datapath,'allD.png'));

totalD = struct('D', Dtot, 'sigD', sigTot);
save(fullfile(datapath,'totalD.mat'), 'totalD');

%% combine gamma

w = 1./(sigGamma.^2);
w = w./sum(w);
gammatot = sum(w.*gamma);
sigTot = 1/sqrt(sum(1./(sigGamma.^2)));

h = errorbar(1:N, gamma, sigGamma,'.');
hold on
errorbar(N+1, gammatot, sigTot,'.r');
hold off
title(['gamma from different movies, gamma1tot = ' num2str(gammatot,'%.2e') '(' num2str(sigTot,'%.2e') ')']);
saveas(gcf, fullfile(datapath,'allGamma.png'));

totalGamma = struct('gamma', gammatot, 'sigGamma', sigTot);
save(fullfile(datapath,'totalGamma.mat'), 'totalGamma');


%% make diffusion constant scatter (old figure draft)

clf
Dvals = [allParam.DperTime(1) 10^(-2) 4*10^(-2) 2*10^(-2) 2*10^(-2)/2];
evals = [allParam.DperTime(1) 2*10^(-3) 10^(-3) 10^(-3) 2*10^(-3)];
lw = 2;
nItems = 5;
cmap = [1 1 0; 0 0 1; 0 1 0; 1 0 0; 0 0 1];
clf
%i = 1;
lw = 4;
marksize = 7;
h = errorbar(1:nItems,Dvals,evals,'o','LineWidth',lw, 'MarkerSize',marksize,'Color','k');
errorbar_tick(h,barwidth,'UNITS');
% lw = 2;
% marksize = 6;
% barwidth = 0.5;
% hold on
% for i = 1:nItems
% h(i) = errorbar(i,Dvals(i),evals(i),'o','LineWidth',lw, 'MarkerSize',marksize,'Color',cmap(i,:));
% errorbar_tick(h(i),barwidth,'UNITS');
% end
% hold off
set(gcf,'color','w');
set(gca,'FontSize', fontsize)
set(gca,'FontWeight', 'bold')
set(gca, 'LineWidth', lw);
set(gca,'XTick', 1:nItems)
set(gca,'XTickLabel',{'Ft-Ds','N-N','Ft','Ds','NCad'})
set(gca,'YTick', 0:0.01:0.2)
ylabel('diffusion constant  (\mum^2/s)');
%title('FRAP measurements of diffusion');
box off
%errorbar(2,10^(-2),10^(-1),'LineWidth',2)
%axis([0 5 0 10^(-2)])