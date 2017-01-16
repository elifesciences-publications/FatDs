%--------------------------------------------------------------------------
%
%   script for analyzing polarity
%
%--------------------------------------------------------------------------

clear all;
close all;

saveResults = true;

% add subfolders of folder containing this script to path
[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

cd(testScriptPath);
cd('..');
addpath(genpath(pwd));

datapath = '/Users/idse/Dropbox/Sprinzak/shared/polarity/fat_ds_rainbow_analysis';
cd(datapath);

files = {...
    fullfile(datapath, '26.8.14 Series070_analysis', '26.8.14 Series070.tif')...
    fullfile(datapath, '26.8.14 Series081_analysis', '26.8.14 Series081.tif')...
    fullfile(datapath, '26.8.14 Series083_analysis', '26.8.14 Series083.tif')...
    fullfile(datapath, '27.1.16_sp8', '008', 'z7', 'Substack (7).tif')...
    fullfile(datapath, '27.1.16_sp8', '008', 'z8', 'Substack (8).tif')...
    fullfile(datapath, '27.1.16_sp8', '008', 'z9', 'Substack (9).tif')...
    fullfile(datapath, '27.1.16_sp8', '14', 'Substack (z5).tif')...
    fullfile(datapath, '27.1.16_sp8', '22', 'z10', 'Substack (10).tif')...
    fullfile(datapath, '27.1.16_sp8', '22', 'z11', 'Substack (11).tif')...
    fullfile(datapath, '27.1.16_sp8', '27', 'z11', 'Substack (11).tif')...
    fullfile(datapath, '27.1.16_sp8', '27', 'z12', 'Substack (12).tif')...
    fullfile(datapath, '27.1.16_sp8', '27', 'z14', 'Substack (14).tif')...
    fullfile(datapath, '27.1.16_sp8', '33', 'z6', 'Substack (6).tif')...
    fullfile(datapath, '27.1.16_sp8', '33', 'z7', 'Substack (7).tif')...
    fullfile(datapath, '27.1.16_sp8', '33', 'z12', 'Substack (12).tif')...
    fullfile(datapath, '27.1.16_sp8', '43', 'z9', 'Substack (9).tif')...
    fullfile(datapath, '27.1.16_sp8', '43', 'z11', 'Substack (11).tif')...
    fullfile(datapath, '27.1.16_sp8', '46', 'z4', 'Substack (4).tif')...
    fullfile(datapath, '27.1.16_sp8', '46', 'z8', 'Substack (8).tif')...
    fullfile(datapath, '27.1.16_sp8', '46', 'z11', 'Substack (11).tif')...
    fullfile(datapath, '27.1.16_sp8', '74', 'z5', 'Substack (5).tif')...
    fullfile(datapath, '27.1.16_sp8', '74', 'z8', 'Substack (8).tif')...
    fullfile(datapath, '27.1.16_sp8', '74', 'z10', 'Substack (10).tif')...
    fullfile(datapath, '1.2.16 f d in the same cell', '17', 'z27', 'Substack (27).tif')...
    fullfile(datapath, '1.2.16 f d in the same cell', '22', 'z21', 'Substack (21).tif')...
    fullfile(datapath, '1.2.16 f d in the same cell', '22', 'z22', 'Substack (22).tif')...
    fullfile(datapath, '1.2.16 f d in the same cell', '22', 'z28', 'Substack (28).tif')...
    fullfile(datapath, '1.2.16 f d in the same cell', '25', 'z25', 'Substack (25).tif')...
    fullfile(datapath, '1.2.16 f d in the same cell', '25', 'z43', 'Substack (43).tif')...
    fullfile(datapath, '1.2.16 f d in the same cell', '30', 'z17', 'Substack (17).tif')...
    fullfile(datapath, '1.2.16 f d in the same cell', '30', 'z20', 'Substack (20).tif')...
    };

datapath = '/Users/idse/Dropbox/Sprinzak/shared/polarity/monoculture rainbows 31.7.16';
files = [files, {...
        fullfile(datapath, '13.7.16 monocul', '17', 'z15','13.7.16 monocul.lif - Series017-z15.tif'),...
        fullfile(datapath, '13.7.16 monocul', '17', 'z18','13.7.16 monocul.lif - Series017-z18.tif'),...
        fullfile(datapath, '13.7.16 monocul', '18', 'z2','17.7.16 monocul.lif - Series018-z2.tif'),...
        fullfile(datapath, '13.7.16 monocul', '20', 'z20','13.7.16 monocul.lif - Series020-z20.tif'),...
        fullfile(datapath, '13.7.16 monocul', '24', 'z17','13.7.16 monocul.lif - Series024-z17.tif'),...
        fullfile(datapath, '13.7.16 monocul', '33', 'z16','13.7.16 monocul.lif - Series033-z16.tif'),...
        fullfile(datapath, '13.7.16 monocul', '35', 'z11','13.7.16 monocul.lif - Series035-z11.tif'),...
        fullfile(datapath, '13.7.16 monocul', '35', 'z17','13.7.16 monocul.lif - Series035-z17.tif'),...
        fullfile(datapath, '13.7.16 monocul', '38', 'z15','13.7.16 monocul.lif - Series038-z15.tif'),...
        fullfile(datapath, '13.7.16 monocul', '41', 'z17','13.7.16 monocul.lif - Series041-z17.tif'),...
        fullfile(datapath, '13.7.16 monocul', '43', 'z15','13.7.16 monocul.lif - Series043-z15.tif'),...
        fullfile(datapath, '13.7.16 monocul', '46', 'z15','13.7.16 monocul.lif - Series046-z15.tif'),...
        fullfile(datapath, '13.7.16 monocul', '50', 'z10','13.7.16 monocul.lif - Series050-z10.tif'),...
        fullfile(datapath, '13.7.16 monocul', '50', 'z25','13.7.16 monocul.lif - Series050-z25.tif'),...
        fullfile(datapath, '13.7.16 monocul', '60', 'z6','13.7.16 monocul.lif - Series060-z6.tif'),...
        fullfile(datapath, '13.7.16 monocul', '84', 'z13','13.7.16 monocul.lif - Series084-z13.tif'),...
        fullfile(datapath, '13.7.16 monocul', '84', 'z23','13.7.16 monocul.lif - Series084-z23.tif'),...
        fullfile(datapath, '13.7.16 monocul', '103', 'z11','13.7.16 monocul.lif - Series103-z11.tif'),...
        fullfile(datapath, '20.7.16 monocul', '7', 'z14','20.7.16 monocul.lif - Series007-z14.tif'),...
        fullfile(datapath, '20.7.16 monocul', '24', 'z3','20.7.16 monocul.lif - Series024-z3.tif'),...
        fullfile(datapath, '20.7.16 monocul', '59', 'z18','20.7.16 monocul.lif - Series059-z18.tif'),...
        fullfile(datapath, '20.7.16 monocul', '70', 'z2','20.7.16 monocul.lif - Series070-z2.tif'),...
        fullfile(datapath, '20.7.16 monocul', '70', 'z4','20.7.16 monocul.lif - Series070-z4.tif'),...
        fullfile(datapath, '20.7.16 monocul', '70', 'z21','20.7.16 monocul.lif - Series070-z21.tif'),...
        fullfile(datapath, '20.7.16 monocul', '80', 'z8','20.7.16 monocul.lif - Series080-z8.tif'),...
        fullfile(datapath, '20.7.16 monocul', '85', 'z10','20.7.16 monocul.lif - Series085-z10.tif')...
        }];
%%
%RGB images are ok here because it is 8-bit data anyway

datapath = '/Users/idse/Dropbox/Sprinzak/shared/polarity/coculture_rainbow_analysis';
cd(datapath);
files = {...
    fullfile(datapath, '7.8.14 coculture series006_analysis', '7.8.14 Series006.tif')...
    fullfile(datapath, '30.6.15 fat and ds co culture','17','30.6.15 fat ds co cul.lif - Series017.tif (RGB).tif')...
    fullfile(datapath, '30.6.15 fat and ds co culture','20','30.6.15 fat ds co cul.lif - Series020.tif (RGB).tif')...
    fullfile(datapath, '30.6.15 fat and ds co culture','29','30.6.15 fat ds co cul.lif - Series029-1.tif')...
    fullfile(datapath, '1.7.15','7','z1','1.7.15 fat ds co cul.lif - Series007-1.tif (RGB).tif'),...
    fullfile(datapath, '1.7.15','7','z5','1.7.15 fat ds co cul.lif - Series007-1.tif (RGB).tif'),...
    fullfile(datapath, '1.7.15','16','z6','1.7.15 fat ds co cul.lif - Series016-1.tif (RGB).tif'),...
    fullfile(datapath, '1.7.15','23','z8','1.7.15 fat ds co cul.lif - Series023-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','15','z5','10.5.15 fat and ds co cul.lif - Series015-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','20','z6','10.5.15 fat and ds co cul.lif - Series020-RGB.tif'),...
    fullfile(datapath, '10.5.16','24','z8','10.5.15 fat and ds co cul.lif - Series024-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','28','z5','10.5.15 fat and ds co cul.lif - Series028-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','36','z7','10.5.15 fat and ds co cul.lif - Series036-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','37','z10','10.5.15 fat and ds co cul.lif - Series037-RGB.tif'),...
    fullfile(datapath, '10.5.16','44','Z3','10.5.15 fat and ds co cul.lif - Series044-RGB.tif'),...
    fullfile(datapath, '10.5.16','44','Z4','10.5.15 fat and ds co cul.lif - Series044-RGB Z4.tif'),...
    fullfile(datapath, '10.5.16','46','z5','10.5.15 fat and ds co cul.lif - Series046-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','47','Z6','10.5.15 fat and ds co cul.lif - Series047-RGB.tif'),...
    fullfile(datapath, '10.5.16','49','z6','10.5.15 fat and ds co cul.lif - Series049-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','50','z5','10.5.15 fat and ds co cul.lif - Series050-RGB.tif'),...
    fullfile(datapath, '10.5.16','54','z7','10.5.15 fat and ds co cul.lif - Series054-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','54','z8','10.5.15 fat and ds co cul.lif - Series054-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','57','z7','10.5.15 fat and ds co cul.lif - Series057-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','82','z12','10.5.15 fat and ds co cul.lif - Series082-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','82','z15','10.5.15 fat and ds co cul.lif - Series082-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','85','z9','10.5.15 fat and ds co cul.lif - Series085-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','102','z13','10.5.15 fat and ds co cul.lif - Series102-1.tif (RGB).tif'),...
    fullfile(datapath, '10.5.16','102','z14','10.5.15 fat and ds co cul.lif - Series102-1.tif (RGB).tif'),...
    fullfile(datapath, '26.4.15 fat ds co cul 60x','13','z10','26.4.15 fat ds co cul 60x.lif - Series063-1.tif (RGB).tif'),...
    fullfile(datapath, '26.4.15 fat ds co cul 60x','24','z7','26.4.15 fat ds co cul 60x.lif - Series024-1.tif (RGB).tif'),...
    fullfile(datapath, '26.4.15 fat ds co cul 60x','28','26.4.15 fat ds co cul 60x.lif - Series028-1.tif')...
    fullfile(datapath, '26.4.15 fat ds co cul 60x','32','26.4.15 fat ds co cul 60x.lif - Series032-1.tif')...
    fullfile(datapath, '26.4.15 fat ds co cul 60x','50','z6','26.4.15 fat ds co cul 60x.lif - Series050-1.tif (RGB).tif'),...
    fullfile(datapath, '26.4.15 fat ds co cul 60x','55','z1','26.4.15 fat ds co cul 60x.lif - Series055-1.tif (RGB).tif'),...
    fullfile(datapath, '26.4.15 fat ds co cul 60x','55','z6','26.4.15 fat ds co cul 60x.lif - Series055-1.tif (RGB).tif'),...
    };

datapath = '/Users/idse/Dropbox/Sprinzak/shared/polarity/co culture 20.7.16';
files = [files, {...
    fullfile(datapath, '6.7.16 co cul','4','z8','6.7.16 co cul.lif - Series004-8.tif'),...
    fullfile(datapath, '6.7.16 co cul','4','z15','6.7.16 co cul.lif - Series004-15.tif'),...
    fullfile(datapath, '6.7.16 co cul','5','z7','6.7.16 co cul.lif - Series005-7.tif'),...
    fullfile(datapath, '6.7.16 co cul','7','z20','6.7.16 co cul.lif - Series007-20.tif'),...
    fullfile(datapath, '6.7.16 co cul','8','z11','6.7.16 co cul.lif - Series008-11.tif'),...
    fullfile(datapath, '6.7.16 co cul','11','z21','6.7.16 co cul.lif - Series011-21.tif'),...
    fullfile(datapath, '6.7.16 co cul','12','z15','6.7.16 co cul.lif - Series012-15.tif'),...
    fullfile(datapath, '6.7.16 co cul','14','z28','6.7.16 co cul.lif - Series014-28.tif'),...
    fullfile(datapath, '6.7.16 co cul','17','z19','6.7.16 co cul.lif - Series017-19.tif'),...
    fullfile(datapath, '6.7.16 co cul','19','z11','6.7.16 co cul.lif - Series019-11.tif'),...
    fullfile(datapath, '6.7.16 co cul','20','z20','6.7.16 co cul.lif - Series020-20.tif'),...
    fullfile(datapath, '6.7.16 co cul','20','z16','6.7.16 co cul.lif - Series020-16.tif'),...
    fullfile(datapath, '6.7.16 co cul','24','z22','6.7.16 co cul.lif - Series024-22.tif'),...
    fullfile(datapath, '6.7.16 co cul','24','z11','6.7.16 co cul.lif - Series024-11.tif')...
    fullfile(datapath, '11.7.16 co cul','9','z16','11.7.16 co cul.lif - Series009-16.tif'),...
    fullfile(datapath, '11.7.16 co cul','13','z21','11.7.16 co cul.lif - Series013-21.tif'),...
    fullfile(datapath, '11.7.16 co cul','13','z28','11.7.16 co cul.lif - Series013-28.tif'),...
    fullfile(datapath, '11.7.16 co cul','16','z20','11.7.16 co cul.lif - Series016-20.tif'),...
    fullfile(datapath, '11.7.16 co cul','16','z22','11.7.16 co cul.lif - Series016-22.tif'),...
    fullfile(datapath, '11.7.16 co cul','23','z18','11.7.16 co cul.lif - Series023-18.tif'),...
    fullfile(datapath, '11.7.16 co cul','26','z15','11.7.16 co cul.lif - Series026-15.tif'),...
    fullfile(datapath, '11.7.16 co cul','30','z16','11.7.16 co cul.lif - Series030-16.tif'),...
    fullfile(datapath, '11.7.16 co cul','34','z22','11.7.16 co cul.lif - Series034-22.tif'),...
    fullfile(datapath, '11.7.16 co cul','37','z19','11.7.16 co cul.lif - Series037-19.tif'),...
    fullfile(datapath, '17.7.16 co cul','6','z7','17.7.16 co cul.lif - Series006-7.tif'),...
    fullfile(datapath, '17.7.16 co cul','6','z10','17.7.16 co cul.lif - Series006-10.tif'),...
    fullfile(datapath, '17.7.16 co cul','7','z6','17.7.16 co cul.lif - Series007-6.tif'),...
    fullfile(datapath, '17.7.16 co cul','9','z2','17.7.16 co cul.lif - Series009-2.tif'),...
    fullfile(datapath, '17.7.16 co cul','9','z7','17.7.16 co cul.lif - Series009-7.tif'),...
    fullfile(datapath, '17.7.16 co cul','11','z1','17.7.16 co cul.lif - Series011-1.tif'),...
    fullfile(datapath, '17.7.16 co cul','11','z12','17.7.16 co cul.lif - Series011-12.tif'),...
    fullfile(datapath, '17.7.16 co cul','12','z3','17.7.16 co cul.lif - Series012-3.tif'),...
    fullfile(datapath, '17.7.16 co cul','15','z4','17.7.16 co cul.lif - Series015-4.tif'),...
    fullfile(datapath, '17.7.16 co cul','19','z2','17.7.16 co cul.lif - Series019-2.tif'),...
    fullfile(datapath, '17.7.16 co cul','20','z2','17.7.16 co cul.lif - Series020-2.tif'),...
    fullfile(datapath, '17.7.16 co cul','20','z11','17.7.16 co cul.lif - Series020-11.tif'),...
    fullfile(datapath, '17.7.16 co cul','25','z2','17.7.16 co cul.lif - Series025-2.tif'),...
    fullfile(datapath, '17.7.16 co cul','28','z8','17.7.16 co cul.lif - Series028-8.tif'),...
    fullfile(datapath, '17.7.16 co cul','31','z2','17.7.16 co cul.lif - Series031-2.tif'),...
    }];

%[fname, filepath] = uigetfile('*.tif');

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
thickness = 41;

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


%% ONLY 6 POINTS IN SNAKE 7.1, LOOK INTO THIS, only 48 in 13.5

for fi = 44%1:numel(files)

%-----------------------------
% read the data
%-----------------------------

[filepath,barefname,ext] = fileparts(files{fi});
resultsdir = filepath;
fname = [barefname ext]; 
[data meta] = readStack(filepath, fname);
data = squeeze(data);

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
figure, imshow(data,[]);
colors = hsv(nSnakes);
hold on
for i = 1:nSnakes
    snake = snakes{1, i};
    %plot(snake(:,1), snake(:,2), 'LineWidth', 2, 'Color', colors(i,:))
    plot(snake(:,1), snake(:,2), '-','LineWidth', 1, 'Color', 'c')
    text(mean(snake(:,1)), mean(snake(:,2)) - 100, num2str(i), 'Color', 'c')
end
hold off

I = getframe(gcf);
imwrite(I.cdata, fullfile(resultsdir, 'overview.png'));
saveas(gcf, fullfile(resultsdir, 'overview.fig'));
close
drawnow

%-----------------------------------
% take the bow out of rainbow
%-----------------------------------

rainbow = {};

for i = 1:nSnakes
    
    snake = snakes{1, i};
    fatBdry = broadImprofile(data(:,:,DsC), snake(:,1), snake(:,2), thickness);
    dsBdry = broadImprofile(data(:,:,FatC), snake(:,1), snake(:,2), thickness);
    rainbow{i} = cat(3, mat2gray(dsBdry), mat2gray(fatBdry), 0*fatBdry);
end

% %% 
% % inspect
% 
% i = 4;
% imshow(rainbow{i},[])

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

% %%
% i =1
% plot(pDs{i}(:,end))
% hold on 
% plot(crap{i},'r')
% hold off

end

%% display some fit

i = 2;
x = 15;
[~crap{i}(x), pDs{i}(x,end)] % quality measure

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
L = 0.11;

DsGradAll = [];
FatGradAll = [];
polAll = [];

cutoffres = 80;
 
for fi = 44%1:numel(files)

    [filepath,barefname,ext] = fileparts(files{fi});
    fname = fullfile(filepath, [barefname '_results.mat']);
    S = load(fname);
    S = S.results;

    x = 30;
    disp('-');

    for i = 1:S.nSnakes
        
        if S.resolution < cutoffres

            idx = S.pDs{i}(:,end) > 0 & ~S.crap{i};

            % 5th parameter is C in A e^(-(y-y0)^2/2s^2) + B + C Theta(y-y0)
            % so is sign of gradient
            FatSgn = 1;%sign(mean(S.pFat{i}(idx,5)));

            DsGrad = FatSgn*mean(S.pDs{i}(idx,5))/mean(S.pDs{i}(idx,1));
            FatGrad = FatSgn*mean(S.pFat{i}(idx,5))/mean(S.pFat{i}(idx,1));

            idx(1:5) = 0;
            idx(end-5:end) = 0;

            %idx = S.pDs{i}(:,end) > 0;
            %idx = ~S.crap{i};

            pol = FatSgn*mean(S.pDs{i}(idx,2) - S.pFat{i}(idx,2))*S.resolution;
            %polcol = max(min((pol + polmax)/(2*polmax),1),0);
            %polcol = max(min((abs(pol)-polmin)/(polmax-polmin),1),0);
            %c = colors(round(1 + polcol*511),:);

            t = polmin;
            if pol > t
                sty = '>';
                c = [1 0 0];
            elseif pol < -t
                sty = '<';
                c = [0 0 1];
            else
                sty = 'o';
                c = [0 1 0];
            end
            %scatter(DsGrad, FatGrad,sty,'MarkerEdgeColor',c)
            text(DsGrad + 0.005, FatGrad, [num2str(fi) '.' num2str(i)])

            DsGradAll = [DsGradAll DsGrad];
            FatGradAll = [FatGradAll FatGrad];
            polAll = [polAll pol];
        end
    end
end

L = 0.25;
plot(L*(-1:1),0*(-1:1),'-b');
plot(0*(-1:1),L*(-1:1),'-b');
hold off
axis equal
axis([-L L -L L]);
xlabel('dDs');
ylabel('dFat');

if saveResults 
    fname = fullfile(datapath, ['gradientVsPolarity_res' num2str(cutoffres) '.png']);
    saveas(gcf, fname);
    fname = fullfile(datapath, ['gradientVsPolarity_res' num2str(cutoffres) '.fig']);
    saveas(gcf, fname);
end

%% simplified plot, binned by FatGrad only

aligned = (DsGradAll > 0 & polAll > polmin) | (DsGradAll < 0 & polAll < -polmin);
notAligned = (DsGradAll > 0 & polAll < -polmin) | (DsGradAll < 0 & polAll > polmin);
notPolarized = polAll < polmin & polAll > -polmin;

binEdge = 0:0.01:0.2;
n = histc(FatGradAll,binEdge);
nAligned = histc(FatGradAll(aligned),binEdge);
nNotAligned = histc(FatGradAll(notAligned),binEdge);
nNotPolarized = histc(FatGradAll(notPolarized),binEdge);

% normalized

bar(binEdge + 0.005, nAligned./n)
axis([binEdge(1) binEdge(end) 0 1])
xlabel('dFat')
ylabel('Fraction polarized along dDs');

if saveResults
    saveas(gcf, fullfile(datapath,['alignedFractionHistogram_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedFractionHistogram_res' num2str(cutoffres) '.fig']));
end

% not normalized

bar(binEdge + 0.005, cat(1,nAligned,nNotAligned,nNotPolarized)', 'stacked')
axis([binEdge(1) binEdge(end) 0 14])
xlabel('dFat')
ylabel('number of boundaries');
legend({'Ds aligned','Ds opposed','not polarized'})

if saveResults
    saveas(gcf, fullfile(datapath,['alignedNumberHistogram_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedNumberHistogram_res' num2str(cutoffres) '.fig']));
end

%% simplified plot, binned by FatGrad only per quadrant 

compatible = DsGradAll < 0;

% COMPATIBLE GRADIENTS
%---------------------

n = histc(FatGradAll(compatible),binEdge);
nAligned = histc(FatGradAll(aligned & compatible),binEdge);
nNotAligned = histc(FatGradAll(notAligned & compatible),binEdge);
nNotPolarized = histc(FatGradAll(notPolarized & compatible),binEdge);

% normalized

bar(binEdge + 0.005, nAligned./n)
axis([binEdge(1) binEdge(end) 0 1])
xlabel('dFat')
ylabel('Fraction polarized along dDs');

if saveResults
    saveas(gcf, fullfile(datapath,['alignedFractionHistogramCompatible_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedFractionHistogramCompatible_res' num2str(cutoffres) '.fig']));
end

% not normalized

bar(binEdge + 0.005, cat(1,nAligned,nNotAligned,nNotPolarized)', 'stacked')
axis([binEdge(1) binEdge(end) 0 14])
xlabel('dFat')
ylabel('number of boundaries');
legend({'Ds aligned','Ds opposed','not polarized'})

if saveResults
    saveas(gcf, fullfile(datapath,['alignedNumberHistogramCompatible_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedNumberHistogramCompatible_res' num2str(cutoffres) '.fig']));
end

% INCOMPATIBLE GRADIENTS
%---------------------

n = histc(FatGradAll(~compatible),binEdge);
nAligned = histc(FatGradAll(aligned & ~compatible),binEdge);
nNotAligned = histc(FatGradAll(notAligned & ~compatible),binEdge);
nNotPolarized = histc(FatGradAll(notPolarized & ~compatible),binEdge);

% normalized

bar(binEdge + 0.005, nAligned./n)
axis([binEdge(1) binEdge(end) 0 1])
xlabel('dFat')
ylabel('Fraction polarized along dDs');

if saveResults
    saveas(gcf, fullfile(datapath,['alignedFractionHistogramIncompatible_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedFractionHistogramIncompatible_res' num2str(cutoffres) '.fig']));
end

% not normalized

bar(binEdge + 0.005, cat(1,nAligned,nNotAligned,nNotPolarized)', 'stacked')
axis([binEdge(1) binEdge(end) 0 14])
xlabel('dFat')
ylabel('number of boundaries');
legend({'Ds aligned','Ds opposed','not polarized'})

if saveResults
    saveas(gcf, fullfile(datapath,['alignedNumberHistogramIncompatible_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedNumberHistogramIncompatible_res' num2str(cutoffres) '.fig']));
end

%% simplest plot, alltogether per quadrant 

compatible = DsGradAll < 0;

% COMPATIBLE GRADIENTS
%---------------------

nC = sum(compatible);
nAlignedC = sum(aligned & compatible);
nNotAlignedC = sum(notAligned & compatible);
nNotPolarizedC = sum(notPolarized & compatible);

nI = sum(~compatible);
nAlignedI = sum(aligned & ~compatible);
nNotAlignedI = sum(notAligned & ~compatible);
nNotPolarizedI = sum(notPolarized & ~compatible);

% normalized

bar([nAlignedC./nC nAlignedI./nI])
axis([0.5 2.5 0 1])
ylabel('Fraction polarized along Ds gradient');
xlabel('Gradients');
set(gca,'XTickLabel',{'compatible','incompatible'})

if saveResults
    saveas(gcf, fullfile(datapath,['alignedFractionCollapsed_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedFractionCollapsed_res' num2str(cutoffres) '.fig']));
end

% not normalized

bar(cat(1,[nAlignedC nAlignedI],[nNotAlignedC nNotAlignedI],[nNotPolarizedC nNotPolarizedI])', 'stacked');
%axis([0.5 3 0 14])
ylabel('Number of boundaries');
xlabel('Gradients');
set(gca,'XTickLabel',{'compatible','incompatible'})
legend({'Ds aligned','Ds opposed','not polarized'},'Location','NorthWest')

if saveResults
    saveas(gcf, fullfile(datapath,['alignedNumberCollapsed_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedNumberCollapsed_res' num2str(cutoffres) '.fig']));
end

%%
%-------------------------------------
% visualize the results per boundary
%-------------------------------------

nrows = 3;
ncols = 2;
scrsz = get(0,'ScreenSize');

for fi = 44%1:numel(files)

    [filepath,barefname,ext] = fileparts(files{fi});
    resultsdir = filepath;
    fname = fullfile(filepath, [barefname '_results.mat']);
    S = load(fname);
    results = S.results;

for i = 1:results.nSnakes
    
    crap = results.crap;
    yskew = results.yskew;
    res = results.resolution;
    pFat = results.pFat;
    pDs = results.pDs;
    rainbow = results.rainbow;
    
    if size(rainbow{i},2) > 10
        
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

    if saveResults
        I = getframe(gcf);
        imwrite(I.cdata, fullfile(resultsdir, ['stats' num2str(i) '.png']));
        saveas(gcf, fullfile(resultsdir, ['stats' num2str(i) '.fig']));
        close;
    end
    
    else
        disp(['skipping snake ' num2str(fi) '.' num2str(i)])
    end
end

end
%pFat{i}(x,:)


%%
%-------------------------------------
% combine boundary fit statistics
%-------------------------------------

nrows = 3;
ncols = 2;
scrsz = get(0,'ScreenSize');

GaussWidthAll = {};
GaussWidthMean = {};
SeparationAll = {};
bdryIdx = 0;

Nbdries = 0;

for fi = 1:numel(files)

    [filepath,barefname,ext] = fileparts(files{fi});
    resultsdir = filepath;
    fname = fullfile(filepath, [barefname '_results.mat']);
    S = load(fname);
    results = S.results;

for i = 1:results.nSnakes
    
    if S.results.resolution < cutoffres
        
        Nbdries = Nbdries + 1;
        
        crap = results.crap;
        yskew = results.yskew;
        res = results.resolution;
        pFat = results.pFat;
        pDs = results.pDs;
        rainbow = results.rainbow;

        if size(rainbow{i},2) > 10

            bdryIdx = bdryIdx + 1;
            GaussWidthAll{bdryIdx, 1} = pDs{i}(~crap{i},3)*res;
            GaussWidthAll{bdryIdx, 2} = pFat{i}(~crap{i}, 3)*res;
    %         GaussWidthMean{bdryIdx, 1} = [mean(GaussWidthAll{bdryIdx, DsC})...
    %                                         std(GaussWidthAll{bdryIdx, DsC})];
    %         GaussWidthMean{bdryIdx, 2} = [mean(GaussWidthAll{bdryIdx, FatC})...
    %                                         std(GaussWidthAll{bdryIdx, FatC})];

            SeparationAll{bdryIdx} = abs((pDs{i}(~crap{i},2) - pFat{i}(~crap{i}, 2))*res);
        else
            disp(['skipping snake ' num2str(fi) '.' num2str(i)])
        end
    end
end

end

%% Gauss width combined

DsWidth = cat(1,GaussWidthAll{:,DsC});
FatWidth = cat(1,GaussWidthAll{:,FatC});

legend = {...
        ['Ds ' num2str(mean(DsWidth), 3) '(' num2str(std(DsWidth), 2) ')'],...
        ['Fat ' num2str(mean(FatWidth), 3) '(' num2str(std(FatWidth), 2) ')']};
    
nhist({DsWidth,FatWidth}, 'noerror', 'color', 'lines', 'linewidth', 2, 'titles', legend, 'maxbins',50);
title(['Boundary Gaussian width, N = ' num2str(Nbdries)]);
xlabel('width (nm)');

if saveResults
    I = getframe(gcf);
    imwrite(I.cdata, fullfile(datapath, 'GaussWidthCombined.png'));
    saveas(gcf, fullfile(datapath, 'GaussWidthCombined.fig'));
    %close;
end

%% separation combined

% looks crappy for monoculture but should be done for coculture

separation = cat(1,SeparationAll{:});

nhist({separation}, 'noerror', 'linewidth', 2, 'maxbins',100);
title(['Fat-Ds separation ' num2str(mean(separation), 3) '(' num2str(std(separation),2) '), Nbdries=' num2str(Nbdries)]);
xlabel('distance (nm)');

if saveResults
    I = getframe(gcf);
    imwrite(I.cdata, fullfile(datapath, 'separationCombined.png'));
    saveas(gcf, fullfile(datapath, 'separationCombined.fig'));
    %close;
end