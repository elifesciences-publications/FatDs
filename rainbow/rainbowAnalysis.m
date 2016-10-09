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

datapath = '/Users/idse/Dropbox/Sprinzak/shared/polarity';
cd(datapath);

fatdsdir = 'fat_ds_rainbow_analysis';

files = {...
    fullfile(datapath, '7.8.14 coculture series006_analysis', '7.8.14 Series006.tif')...
    fullfile(datapath, '26.8.14 Series070_analysis', '26.8.14 Series070.tif')...
    fullfile(datapath, '26.8.14 Series081_analysis', '26.8.14 Series081.tif')...
    fullfile(datapath, '26.8.14 Series083_analysis', '26.8.14 Series083.tif')...
    fullfile(datapath, fatdsdir, '27.1.16_sp8', '008', 'z7', 'Substack (7).tif')...
    fullfile(datapath, fatdsdir, '27.1.16_sp8', '008', 'z8', 'Substack (8).tif')...
    fullfile(datapath, fatdsdir, '27.1.16_sp8', '008', 'z9', 'Substack (9).tif')...
    };
%[fname, filepath] = uigetfile('*.tif');

fi = 1;

[filepath,barefname,ext] = fileparts(files{fi});
fname = [barefname ext];
snakeFile = fullfile(filepath, [barefname ' snakes.txt']);

% % make a directory for saving the results
% resultsdir = fullfile(filepath, [barefname '_analysis']);
% if ~exist(resultsdir, 'dir')
%     mkdir(resultsdir);
% end

resultsdir = filepath;

% SETTINGS! 
%-----------------

% channel indices
DsC = 1;
FatC = 2;

% resolution, nm per pixel
res = 60;

% emissions peak wavelengths, 530, 610

%-----------------------------
% read the data
%-----------------------------

[data nChannels] = readStack(filepath, fname);
data = squeeze(data);

%% 
figure,
imshow(data,[])

%%
%-----------------------------
% read snake
%-----------------------------

%snakeFile = uigetfile('*'); % ['snakes' barefname]
%snakeFile = fullfile(filepath, snakeFile);
snakes = readSnake(snakeFile, 1);
nSnakes = size(snakes,2);

% visualize and save
figure, imshow(data,[]);
colors = hsv(nSnakes);
hold on
for i = 1:nSnakes
    snake = snakes{1, i};
    %plot(snake(:,1), snake(:,2), 'LineWidth', 2, 'Color', colors(i,:))
    plot(snake(:,1), snake(:,2), 'LineWidth', 2, 'Color', 'c')
    text(mean(snake(:,1)), mean(snake(:,2)) - 100, num2str(i), 'Color', 'c')
end
hold off

I = getframe(gcf);
imwrite(I.cdata, fullfile(resultsdir, 'overview.png'));
saveas(gcf, fullfile(resultsdir, 'overview.fig'));
%close

%%
%-----------------------------------
% take the bow out of rainbow
%-----------------------------------

rainbow = {};

% SETTINGS
thickness = 51;

for i = 1:nSnakes
    
    snake = snakes{1, i};
    fatBdry = broadImprofile(data(:,:,2), snake(:,1), snake(:,2), thickness);
    dsBdry = broadImprofile(data(:,:,1), snake(:,1), snake(:,2), thickness);
    rainbow{i} = cat(3, mat2gray(dsBdry), mat2gray(fatBdry), 0*fatBdry);
end

%% 
% inspect

i = 1;
imshow(rainbow{i},[])

%%
%----------------------------------
% fit Gaussian
%----------------------------------

% function f to be fitted: A e^(-(y-y0)^2/2s^2) + B + C Theta(y-y0)
% parameters p = [A y0 s B]
f = @(p,y) p(1)*exp(- (y-p(2)).^2/(2*p(3)^2) ) + p(4) + p(5)*(y>p(2));
y = 1:thickness;

% initial estimates for parameters
pinit = [1, thickness/2, 500/(2*res), 0.1, 0];
nParam = length(pinit);

pFat = {};
pDs = {};
crap = {};

for i = 1:nSnakes
    
    % last parameter entry: residuals
    pFat{i} = zeros([thickness nParam + 1]);
    pDs{i} = zeros([thickness nParam + 1]);

    for x = 1:size(rainbow{i},2)

        Ifat = rainbow{i}(:,x,2);
        Ids = rainbow{i}(:,x,1);

        % get the Fat parameters
        % E: energy functional
        E = @(p) sum((Ifat' - f(p,y)).^2);
        [p, fminres] = fminsearch(E, pinit);
        pFat{i}(x,:) = [p, fminres];

        % get the Ds parameters
        E = @(p) sum((Ids' - f(p,y)).^2);
        [p, fminres] = fminsearch(E, pinit);
        pDs{i}(x,:) = [p, fminres];
    end

    % crap mask
    %-----------

    % width unrealistically large
    crap{i} = pFat{i}(:,3) > 10 | pDs{i}(:,3) > 10;
    % amplitude too small
    crap{i} = crap{i} | (pFat{i}(:,1) < 0.1 | pDs{i}(:,1) < 0.1);
    % background too high
    crap{i} = crap{i} | (pFat{i}(:,4) > 0.4 | pDs{i}(:,4) > 0.4);
end

%% display some fit

i = 1;
x = 10;
crap{i}(x) % check if this fit is crap

y = 1:thickness;
plot(squeeze(rainbow{i}(:,x,1)), '-r', 'LineWidth', 1)

hold on;
plot(squeeze(rainbow{i}(:,x,2)), '-g', 'LineWidth', 1)
plot(f(pDs{i}(x,:),y), '-r', 'LineWidth', 2)
plot(f(pFat{i}(x,:),y), '-g', 'LineWidth', 2)

plot(pFat{i}(x,4) + pFat{i}(x,5)*(y>pFat{i}(x,2)),'-b')
plot(pDs{i}(x,4) + pDs{i}(x,5)*(y>pDs{i}(x,2)),'-b')

hold off;

I = getframe(gcf);
imwrite(I.cdata, fullfile(resultsdir, 'someFit.png'));
saveas(gcf, fullfile(resultsdir, 'someFit.fig'));
%close

%% 
x = 30;
disp('-');
L = 0.15;
clf
hold on
for i = 1:nSnakes
    pol = mean(pDs{i}(~crap{i},2) - pFat{i}(~crap{i},2));
    if pol > 0
        sty = '>';
    else
        sty = '<';
    end
    scatter(mean(pDs{i}(~crap{i},5)), mean(pFat{i}(~crap{i},5)),sty)
end
plot(L*(-1:1),0*(-1:1),'-b');
plot(0*(-1:1),L*(-1:1),'-b');
hold off
axis equal
axis([-L L -L L]);
xlabel('dDs');
ylabel('dFat');

%%
%--------------------------------
% measure skewness of profiles 
%--------------------------------

yskew = {};

for c = 1:2
    for i = 1:nSnakes

        x = ~crap{i};
        profilesDs = squeeze(rainbow{i}(:,x, c));
        profilesDs = profilesDs./repmat(sum(profilesDs, 1), [thickness 1]);

        nx = size(profilesDs, 2);
        yavg = zeros([nx 1]);
        ystd = zeros([nx 1]);
        yskew{i, c} = zeros([nx 1]);

        for x = 1:size(profilesDs, 2)

            yavg(x) = sum(y.*profilesDs(:,x)');
            ystd(x) = sqrt(sum((y-yavg(x)).^2.*profilesDs(:,x)'));
            yskew{i, c}(x) = sum((y-yavg(x)).^3.*profilesDs(:,x)')/ystd(x)^3;

        end
    end
end

hist(yskew{1, DsC})


%%
%-----------------------
% visualize the results
%-----------------------

nrows = 3;
ncols = 2;
scrsz = get(0,'ScreenSize');

for i = 1:nSnakes
    
    figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)])
    subplot(nrows,ncols,1)
    bla = rainbow{i};
    bla(:,crap{i},:) = 0.5*bla(:,crap{i},:);
    imshow(bla,[])
    title('Boundary');

    % width statistic
    subplot(nrows,ncols,2)
    data = {pDs{i}(~crap{i},3)*res, pFat{i}(~crap{i},3)*res};
    legend = {...
        ['Ds ' num2str(mean(data{DsC}), 3) '(' num2str(std(data{DsC}), 2) ')'],...
        ['Fat ' num2str(mean(data{FatC}), 3) '(' num2str(std(data{FatC}), 2) ')']};
    nhist(data, 'noerror', 'color', 'lines', 'linewidth', 2, 'titles', legend);
    title('Boundary Gaussian width');
    xlabel('width (nm)');

    % overlay the profiles
    subplot(nrows,ncols,3)
    x = ~crap{i};
    y = (1:thickness)*res;
    plot(y, squeeze(rainbow{i}(:,x,1)), '-r')
    hold on;
    plot(y, squeeze(rainbow{i}(:,x,2)), '-g')
    hold off;
    title('Intensity profiles');
    xlabel('distance (nm)');

    % distance statistics
    subplot(nrows,ncols,4)
    distance = (pDs{i}(~crap{i},2) - pFat{i}(~crap{i}, 2))*res;
    nhist({distance}, 'noerror')
    mean(distance)
    std(distance)
    title(['Fat-Ds separation ' num2str(mean(distance), 3) '(' num2str(std(distance),2) ')']);
    xlabel('distance (nm)');
    
    subplot(nrows,ncols,5)
    legend = {'Ds', 'Fat'};
    nhist(yskew(i,:), 'noerror', 'color', 'lines', 'linewidth', 2, 'titles', legend);
    title('Profile skewness');

    I = getframe(gcf);
    imwrite(I.cdata, fullfile(resultsdir, ['stats' num2str(i) '.png']));
    saveas(gcf, fullfile(resultsdir, ['stats' num2str(i) '.fig']));
end

%pFat{i}(x,:)




%%

% FIT DOUBLE PEAK IN ONE CHANNEL?
