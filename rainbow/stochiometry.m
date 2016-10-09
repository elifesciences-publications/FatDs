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
    };

zidx = {1,3,2,1,3,3,2,1,1,2};

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

for fi = 1%:numel(files)

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

%% get boundary intensities from each snake

bdryFatTot = [];
bdryDsTot = [];
bdryFatMean = [];
bdryDsMean = [];
fatBdry = {};
dsBdry = {};

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

    % %% 
    % figure,
    % imshow(data,[])

    %-----------------------------
    % read snake
    %-----------------------------

    %snakeFile = uigetfile('*'); % ['snakes' barefname]
    %snakeFile = fullfile(filepath, snakeFile);

    snakeFile = fullfile(filepath, [barefname ' snakes.txt']);
    if ~exist(snakeFile,'file')
        snakeFile = fullfile(filepath, 'snakes.txt');
    end

    snakes = readSnake(snakeFile, 1);
    nSnakes = size(snakes,2);
    %nSnakes = 2;

    % visualize and save
    figure, imshow(cat(3,imadjust(data(:,:,1)),imadjust(data(:,:,2)),0*data(:,:,3)),[]);
    colors = hsv(nSnakes);
    hold on
    for i = 1:nSnakes
        snake = snakes{1, i};
        %plot(snake(:,1), snake(:,2), 'LineWidth', 2, 'Color', colors(i,:))
        plot(snake(:,1), snake(:,2), '-','LineWidth', 1, 'Color', 'c')
        text(mean(snake(:,1)), mean(snake(:,2)) - 100, num2str(i), 'Color', 'c')
    end
    hold off

    %I = getframe(gcf);
    %imwrite(I.cdata, fullfile(resultsdir, 'overview.png'));
    %saveas(gcf, fullfile(resultsdir, 'overview.fig'));
    close
    drawnow

    %-----------------------------------
    % take the bow out of rainbow
    %-----------------------------------

    rainbow = {};

    for i = 1:nSnakes


        snake = snakes{1, i};
        fatBdry = [fatBdry, broadImprofile(data(:,:,2), snake(:,1), snake(:,2), thickness)];
        dsBdry = [dsBdry, broadImprofile(data(:,:,1), snake(:,1), snake(:,2), thickness)];

        %minI = min(min(dsBdry{end}(:)), min(fatBdry{end}(:)));
        %maxI = max(max(dsBdry{end}(:)), max(fatBdry{end}(:)));
        
        % rescaling rainbow LUT happens here 
        %---------------------------------------
        % fixed range so parameters have meaning
        minI = 200;
        maxI = 10^4;
        rainbow{i} = cat(3, mat2gray(dsBdry{end},[minI maxI]), mat2gray(fatBdry{end},[minI maxI]), 0*fatBdry{end});

        bdryFatTot = [bdryFatTot sum(fatBdry{end}(:))];
        bdryDsTot = [bdryDsTot sum(dsBdry{end}(:))];
        
        bdryFatMean = [bdryFatMean mean(fatBdry{end}(:))];
        bdryDsMean = [bdryDsMean mean(dsBdry{end}(:))];
    end


    %----------------------------------
    % fit Gaussian
    %----------------------------------

    % function f to be fitted: A e^(-(y-y0)^2/2s^2) + B + C Theta(y-y0)
    % parameters p = [A y0 s B C]
    f = @(p,y) p(1)*exp(- (y-p(2)).^2/(2*p(3)^2) ) + p(4) + p(5)*(y>p(2));
    y = 1:thickness;

    % initial estimates for parameters
    pinit = [1, thickness/2, 500/(2*res), 0.1, 0];
    nParam = length(pinit);

    pFat = {};
    pDs = {};
    crap = {};
    yskew = {};

    for i = 1:nSnakes

        snakeLength = size(snakes{1,i}, 1);

        % last 2 parameter entries: residuals, exitflag
        pFat{i} = zeros([snakeLength nParam + 2]);
        pDs{i} = zeros([snakeLength nParam + 2]);

        for x = 1:size(rainbow{i},2)

            Ifat = rainbow{i}(:,x,2);
            Ids = rainbow{i}(:,x,1);

            % get the Fat parameters
            % E: energy functional
            E = @(p) sum((Ifat' - f(p,y)).^2);
            [p, fminres, exitflag] = fminsearch(E, pinit);
            pFat{i}(x,:) = [p, fminres, exitflag];

            % get the Ds parameters
            E = @(p) sum((Ids' - f(p,y)).^2);
            [p, fminres, exitflag] = fminsearch(E, pinit);
            pDs{i}(x,:) = [p, fminres, exitflag];
        end

        % crap mask
        %-----------

        % width unrealistically large
        crap{i} = pFat{i}(:,3)*res > 500 | pDs{i}(:,3)*res > 500;
        % amplitude too small
        crap{i} = crap{i} | (pFat{i}(:,1) < 0.1 | pDs{i}(:,1) < 0.1);
        % background too high
        crap{i} = crap{i} | (pFat{i}(:,4) > 0.4 | pDs{i}(:,4) > 0.4);
        % bad fit (e.g. vesicle nearby)
        crap{i} = crap{i} | pDs{i}(:,6) > 0.5 | pFat{i}(:,6) > 0.5;

        %--------------------------------
        % measure skewness of profiles 
        %--------------------------------

        for c = 1:2

            x = ~crap{i};
            profiles = squeeze(rainbow{i}(:,x, c));
            profiles = profiles./repmat(sum(profiles, 1), [thickness 1]);

            nx = size(profiles, 2);
            yavg = zeros([nx 1]);
            ystd = zeros([nx 1]);
            yskew{i, c} = zeros([nx 1]);

            for x = 1:size(profiles, 2)

                yavg(x) = sum(y.*profiles(:,x)');
                ystd(x) = sqrt(sum((y-yavg(x)).^2.*profiles(:,x)'));
                yskew{i, c}(x) = sum((y-yavg(x)).^3.*profiles(:,x)')/ystd(x)^3;

            end
        end

        %hist(yskew{1, DsC})
    end

    % save results

    % function f to be fitted: A e^(-(y-y0)^2/2s^2) + B + C Theta(y-y0)
    % parameters p = [A y0 s B C]
    results = struct('nSnakes', nSnakes, 'pDs',{pDs}, 'pFat',{pFat},...
                        'crap', {crap}, 'resolution', res, 'yskew',{yskew},...
                        'rainbow', {rainbow});
    fname = fullfile(filepath, [barefname '_results.mat']);
    save(fname, 'results')

end

%% 
%---------------------------------------
% combine all results in scatter plot
%---------------------------------------

% function f to be fitted: A e^(-(y-y0)^2/2s^2) + B + C Theta(y-y0)
% parameters p = [A y0 s B C]

% ERROR BARS
figure,
clf 
hold on 
colors = jet(512);
polmax = 200;
polmin = 20;
L = 1;

DsAmpAll = [];
FatAmpAll = [];
label = {};

for fi = 1:numel(files)

    [filepath,barefname,ext] = fileparts(files{fi});
    fname = fullfile(filepath, [barefname '_results.mat']);
    S = load(fname);
    S = S.results;

    x = 30;
    disp('-');

    for i = 1:S.nSnakes
        
        %if S.resolution < 80
                
        idx = S.pDs{i}(:,end) > 0 & ~S.crap{i};
        
        FatSgn = sign(mean(S.pFat{i}(idx,5)));
        
        DsAmp = median(S.pDs{i}(idx,1));
        FatAmp = median(S.pFat{i}(idx,1));

        DsAmpAll = [DsAmpAll DsAmp];
        FatAmpAll = [FatAmpAll FatAmp];
        
        idx(1:5) = 0;
        idx(end-5:end) = 0;
        
        %idx = S.pDs{i}(:,end) > 0;
        %idx = ~S.crap{i};
        
        pol = FatSgn*mean(S.pDs{i}(idx,2) - S.pFat{i}(idx,2))*S.resolution;
        polcol = max(min((pol + polmax)/(2*polmax),1),0);
        %polcol = max(min((abs(pol)-polmin)/(polmax-polmin),1),0);
        c = colors(round(1 + polcol*511),:);

        scatter(DsAmp, FatAmp,'.','MarkerEdgeColor','r')
        text(DsAmp + 0.005, FatAmp, [num2str(fi) '.' num2str(i)])
        
        label = [label [num2str(fi) '.' num2str(i)]];
        %end
    end
end

hold off
axis equal
%axis([0 L 0 L]);
xlabel('Ds amplitude');
ylabel('Fat amplitude');

C = corrcoef(FatAmpAll,DsAmpAll);
title(['corr coeff ' num2str(C(1,2),2)]);

% amplitudes are scaled with fixed lookup table so 
% minI + A*(maxI - minI) is true amplitude
% the reason the amplitudes are bad is that they are 230 nm resolution

fname = fullfile(datapath, 'FatAmpVsDsAmpMedian.png');
saveas(gcf, fname);

save(fullfile(datapath, 'fitAmplitudes'),'FatAmpAll','DsAmpAll');

%% 
%---------------------------------------
% make scatter plot
%---------------------------------------
figure,
scatter(bdryFatTot, bdryDsTot,'.');
idx = 1:numel(bdryFatTot);
for i = idx 
    text(bdryFatTot(i)+10^5, bdryDsTot(i), label(i));
end
xlabel('Fat intensity');
ylabel('Ds intensity');
axis equal
axis([0 max(bdryFatTot(:)) 0 max(bdryFatTot(:))]);

C = corrcoef(bdryFatTot,bdryDsTot);
title(['Total intensity snake corr coeff ' num2str(C(1,2),2)]);

saveas(gcf, fullfile(datapath, 'FatVsDsTotSnake.png'));
save(fullfile(datapath, 'totalBdryI'),'bdryFatTot','bdryDsTot','bdryFatMean','bdryDsMean');

% should not use total!
% of course a large boundary has both more total green and more total red,
% that correlation is simply an effect of size

%% correlate means

figure,
scatter(bdryFatMean, bdryDsMean,'.');
idx = 1:numel(bdryFatMean);
for i = idx 
    text(bdryFatMean(i), bdryDsMean(i), label(i));
end
xlabel('Fat intensity');
ylabel('Ds intensity');
axis equal
axis([0 max(bdryFatMean(:)) 0 max(bdryFatMean(:))]);

C = corrcoef(bdryFatMean,bdryDsMean);
title(['Total intensity snake corr coeff ' num2str(C(1,2),2)]);

saveas(gcf, fullfile(datapath, 'FatVsDsMeanSnake.png'));

% SEE MORE IN MOVIE ANALYSIS
