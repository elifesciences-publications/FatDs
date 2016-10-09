%% plot of values from different movies

FRAPdataLocationsNCad
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

