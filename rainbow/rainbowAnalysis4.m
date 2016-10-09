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

% SETTINGS! 
%-----------------

% % resolution, nm per pixel
% res = 60;

% SETTINGS
thickness = 31;

% emissions peak wavelengths, 530, 610

%MONOCULTURE
%channel indices
DsC = 1;
FatC = 2;
DAPIC = 3;
 
oldmonofiles = {...
    fullfile(datapath, '27.1.16_sp8', '008', 'z7', 'Substack (7).tif')...
    fullfile(datapath, '27.1.16_sp8', '008', 'z8', 'Substack (8).tif')...
    fullfile(datapath, '27.1.16_sp8', '008', 'z9', 'Substack (9).tif')...
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

%%
DsC = 2;
FatC = 1;
DAPIC = 3;
datapath = '/Users/idse/Dropbox/Sprinzak/shared/polarity/monoculture rainbows 31.7.16';
files = [oldmonofiles, {...
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
% channel indices
DsC = 2;
FatC = 1;
DAPIC = 3;

% COCULTURE
datapath = '/Users/idse/Dropbox/Sprinzak/shared/polarity/co culture 20.7.16';
files = [{...
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

for fi = 1:numel(files)

%-----------------------------
% read the data
%-----------------------------

[filepath,barefname,ext] = fileparts(files{fi});
resultsdir = filepath;
fname = [barefname ext]; 
[data meta] = readStack(filepath, fname);
data = squeeze(data);

%MONOCULTURE:
if i <= numel(oldmonofiles)
    % for all files in polarity/fat_ds_rainbow_analysis;
    DsC = 1;
    FatC = 2;
else
   % all files in polarity/monoculture rainbows 31.7.16
   DsC = 2;
   FatC = 1;
end
if strfind(fname,'Substack')
    
    s = strsplit(fname,{'(',')'});
    zslice = str2double(s{2});
    s = strsplit(files{fi},'/');
    series = str2double(s{end-2});
    fullstackfname = sprintf('27.1.16 fat ds in the same cell.lif - Series%.3d.tif',series);
    if ~exist(fullfile(filepath,'..',fullstackfname),'file')
        fullstackfname = sprintf('1.2.16 fat ds in the same cell.lif - Series%.3d.tif',series);
    end
else
    s = strsplit(fname,'-','CollapseDelimiters',false);
    fullstackfname = s{1};
    for i = 2:numel(s)-1
        fullstackfname = [fullstackfname '-' s{i}];
    end
    fullstackfname = [fullstackfname '.tif'];
    s = strsplit(s{end},'.');
    zslice = str2double(s{1}(2:end));
end

% % COCULTURE:
% % dumb stuff to try MIP, SIP and basal gradient readout
% s = strsplit(fname,'-','CollapseDelimiters',false);
% fullstackfname = s{1};
% for i = 2:numel(s)-1
%     fullstackfname = [fullstackfname '-' s{i}];
% end
% fullstackfname = [fullstackfname '.tif'];
% s = strsplit(s{end},'.');
% zslicestr = s{1};
% zslice = str2double(zslicestr);

stack = readStack(fullfile(filepath,'..'),fullstackfname);
basalz = min(zslice+2,size(stack,3));

DAPIMIP = mean(stack(:,:,zslice:basalz,DAPIC),3); %sum(stack(:,:,:,DAPIC),3);
fatMIP = mean(stack(:,:,zslice:basalz,FatC),3);%sum(stack(:,:,:,FatC),3);
dsMIP = mean(stack(:,:,zslice:basalz,DsC),3);%sum(stack(:,:,:,DsC),3);

res = round(1000*meta.xres); % resolution in nm/pixel

% figure,
% imshofiw(data,[])

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

%-----------------------------------
% take the bow out of rainbow
%-----------------------------------

rainbow = {};
dsCyt = {};
fatCyt = {};
nx = {};
ny = {};

for i = 1:nSnakes
    
    % for the rainbow, read out one slice
    d = 30;
    snake = snakes{1, i};
    [fatBdry,~,~,nx{i},ny{i}] = broadImprofile(data(:,:,FatC), snake(:,1), snake(:,2), thickness + 2*d);
    dsBdry = broadImprofile(data(:,:,DsC), snake(:,1), snake(:,2), thickness + 2*d);
    DAPIBdry = broadImprofile(data(:,:,DAPIC), snake(:,1), snake(:,2), thickness + 2*d);
    rainbow{i} = cat(3, mat2gray(dsBdry(d+1:end-d,:)),...
                            mat2gray(fatBdry(d+1:end-d,:)), 0*fatBdry(d+1:end-d,:));

    % FOR 13.7.16, originally d=50, fth = dth = 35, zbasal = zslice + 4
    % now for 13.7.16 and 20.7.16, d=30 fth = dth = 50, zbasal = zslice + 2
                        
    r=5;
%     if i < 29
%        fth = 50;
%        dth = 50;
%     else
        fth = 50;
        dth = 50;
%     end


    % for the polarity, read out several slices
    fatBdry = broadImprofile(fatMIP, snake(:,1), snake(:,2), thickness + 2*d);
    dsBdry = broadImprofile(dsMIP, snake(:,1), snake(:,2), thickness + 2*d);
    DAPIBdry = broadImprofile(DAPIMIP, snake(:,1), snake(:,2), thickness + 2*d);
    
    DAPImask = imfill(imclose(imopen(DAPIBdry > 10,strel('disk',2)),strel('disk',3)),'holes');
    fatMask = imdilate(fatBdry > fth, strel('disk',r));
    dsMask = imdilate(dsBdry > dth, strel('disk',r));
    mask = ~(DAPImask | fatMask | dsMask);
    Amask = mask(1:d,:);
    Bmask = mask(end-d+1:end,:);

    A = fatBdry(1:d,:);
    B = fatBdry(end-d+1:end,:);
    %Amask = imerode(A < fth, strel('disk',r));
    %Bmask = imerode(B < fth, strel('disk',r));
    fatCyt{i} = [mean(A(Amask)) mean(B(Bmask)) median(A(Amask)) median(B(Bmask)) (sum(Amask(:))+sum(Bmask(:)))/(2*numel(Amask))];
    disp([num2str(i) ' F: ' num2str(fatCyt{i}(1)/fatCyt{i}(2),2) ';' num2str(fatCyt{i}(3)/fatCyt{i}(4),2)]);

    A = dsBdry(1:d,:);
    B = dsBdry(end-d+1:end,:);
    %Amask = imerode(A < dth, strel('disk',r));
    %Bmask = imerode(B < dth, strel('disk',r));
    dsCyt{i} = [mean(A(Amask)) mean(B(Bmask)) median(A(Amask)) median(B(Bmask)) (sum(Amask(:))+sum(Bmask(:)))/(2*numel(Amask))];
    disp([num2str(i) ' D: ' num2str(dsCyt{i}(1)/dsCyt{i}(2),2) ';' num2str(dsCyt{i}(3)/dsCyt{i}(4),2)]);
end

%figure, imshow(imadjust(mat2gray(A)).*Amask,[])

% %% 
% % inspect
% 
% i = 4;

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
    
%     %--------------------------------
%     % measure skewness of profiles 
%     %--------------------------------
% 
%     for c = 1:2
% 
%         x = ~crap{i};
%         profiles = squeeze(rainbow{i}(:,x, c));
%         profiles = profiles./repmat(sum(profiles, 1), [thickness 1]);
% 
%         nx = size(profiles, 2);
%         yavg = zeros([nx 1]);
%         ystd = zeros([nx 1]);
%         yskew{i, c} = zeros([nx 1]);
% 
%         for x = 1:size(profiles, 2)
% 
%             yavg(x) = sum(y.*profiles(:,x)');
%             ystd(x) = sqrt(sum((y-yavg(x)).^2.*profiles(:,x)'));
%             yskew{i, c}(x) = sum((y-yavg(x)).^3.*profiles(:,x)')/ystd(x)^3;
% 
%         end
%     end

    %hist(yskew{1, DsC})
end

% save results
%------------------
% function f to be fitted: A e^(-(y-y0)^2/2s^2) + B + C Theta(y-y0)
% parameters p = [A y0 s B C]
results = struct('nSnakes', nSnakes, 'pDs',{pDs}, 'pFat',{pFat},...
                    'crap', {crap}, 'resolution', res, 'yskew',{yskew},...
                    'rainbow', {rainbow},...
                    'fatCyt', {fatCyt}, 'dsCyt', {dsCyt});
fname = fullfile(filepath, [barefname '_results.mat']);
save(fname, 'results')

% visualize and save visualization
%----------------------------------
figure, imshow(cat(3,   imadjust(mat2gray(data(:,:,DsC))),...
                        imadjust(mat2gray(data(:,:,FatC))),mat2gray(data(:,:,3))),[]);
lw = 2;
l = 30;
s = 4;
polmin = 20;
polmax = 200;
hold on
for i = 1:nSnakes
    
    % gradient signs 
    S = results;
    idx = S.pDs{i}(:,end) > 0 & ~S.crap{i};
    j = 1; % THIS DECIDES BETWEEN MEDIAN AND MEAN!!
    FatSgn = 2*(double(S.fatCyt{i}(j+1) > S.fatCyt{i}(j)) - 1/2);
    DsSgn = 2*(double(S.dsCyt{i}(j+1) > S.dsCyt{i}(j)) - 1/2);

    % polarization
    idx = S.pDs{i}(:,end) > 0 & ~S.crap{i};
    pol = -mean(S.pDs{i}(idx,2) - S.pFat{i}(idx,2))*S.resolution;
    polsgn = pol/abs(pol);
    
    snake = snakes{1, i};
    %plot(snake(:,1), snake(:,2), 'LineWidth', 2, 'Color', colors(i,:))
    if abs(pol) > polmin && abs(pol) < polmax
        plot(snake(:,1), snake(:,2), '-','LineWidth', lw, 'Color', 'r')
        plot(snake(:,1) + s*polsgn*nx{i}', snake(:,2) + s*polsgn*ny{i}', '-','LineWidth', lw, 'Color', 'g')
    else
        plot(snake(:,1), snake(:,2), '-','LineWidth', lw, 'Color', 'y')
    end
    
    text(mean(snake(:,1)), mean(snake(:,2)) - 50, num2str(i), 'Color', 'c')
    snakecen = round(size(snake,1)/2);
    quiver(snake(snakecen,1),snake(snakecen,2),l*FatSgn*nx{i}(snakecen),l*FatSgn*ny{i}(snakecen),'g','LineWidth',2);
    quiver(snake(snakecen+5,1),snake(snakecen+5,2),l*DsSgn*nx{i}(snakecen),l*DsSgn*ny{i}(snakecen),'r','LineWidth',2);
end
hold off
saveas(gcf, fullfile(resultsdir, 'overview2.png'));
saveas(gcf, fullfile(resultsdir, 'overview2.fig'));
close
drawnow

% %%
% i =1
% plot(pDs{i}(:,end))
% hold on 
% plot(crap{i},'r')
% hold off

end

%% display some fit
figure,
i = 3;
x = 10;
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

% are gradients of Fat and Ds the same or opposing ?
Ndifferent = 0;
Nsame = 0;
Nall = 0;

for fi = 1:numel(files)

    [filepath,barefname,ext] = fileparts(files{fi});
    fname = fullfile(filepath, [barefname '_results.mat']);
    S = load(fname);
    S = S.results;

    x = 30;
    %disp('-');

    for i = 1:S.nSnakes
        
        Nall = Nall + 1;
        
        idx = S.pDs{i}(:,end) > 0 & ~S.crap{i};
        FatSgn = 2*(double(S.fatCyt{i}(2) > S.fatCyt{i}(1)) - 1/2);
        j = 1; % j=1 -> mean, j=3-> median value of box away bdry
        DsGrad = FatSgn*(S.dsCyt{i}(j+1) - S.dsCyt{i}(j))/(S.dsCyt{i}(j+1) + S.dsCyt{i}(j));
        FatGrad = FatSgn*(S.fatCyt{i}(j+1) - S.fatCyt{i}(j))/(S.fatCyt{i}(j+1) + S.fatCyt{i}(j));

        if S.resolution < cutoffres...% && ~(any(S.dsCyt{i}(1:2) <= 1.5) && abs(DsGrad) < 0.25)...
            && ~(S.dsCyt{i}(5) < 0.1)%|| S.fatCyt{i}(5) < 0.9)
            %&& ~(any(S.dsCyt{i}(3:4) == 0) || any(S.fatCyt{i}(3:4) == 0))
            %&& ~(S.dsCyt{i}(5) < 0.5)%|| S.fatCyt{i}(5) < 0.9)
            
            %disp(S.dsCyt{i});
            %disp(S.fatCyt{i});
            
            if DsGrad > 0
                Nsame = Nsame + 1;
            else
                Ndifferent = Ndifferent + 1;
            end
            
            %idx(1:5) = 0;
            %idx(end-5:end) = 0;

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
            scatter(DsGrad, FatGrad,sty,'MarkerEdgeColor',c)
            %text(DsGrad + 0.005, FatGrad, [num2str(fi) '.' num2str(i)])

            DsGradAll = [DsGradAll DsGrad];
            FatGradAll = [FatGradAll FatGrad];
            polAll = [polAll pol];
        end
    end
end

L = 1;
plot(L*(-1:1),0*(-1:1),'-b');
plot(0*(-1:1),L*(-1:1),'-b');
hold off
axis equal
axis([-L L 0 L]);
xlabel('dDs');
ylabel('dFat');

if saveResults 
    fname = fullfile(datapath, ['gradientVsPolarityNew2_res' num2str(cutoffres) '.png']);
    saveas(gcf, fname);
    fname = fullfile(datapath, ['gradientVsPolarityNew2_res' num2str(cutoffres) '.fig']);
    saveas(gcf, fname);
end

Nsame
Ndifferent

%% simplified plot, binned by FatGrad only

aligned = (DsGradAll > 0 & polAll > polmin) | (DsGradAll < 0 & polAll < -polmin);
notAligned = (DsGradAll > 0 & polAll < -polmin) | (DsGradAll < 0 & polAll > polmin);
notPolarized = polAll < polmin & polAll > -polmin;

binEdge = 0:0.1:0.7;
n = histc(FatGradAll,binEdge);
nAligned = histc(FatGradAll(aligned),binEdge);
nNotAligned = histc(FatGradAll(notAligned),binEdge);
nNotPolarized = histc(FatGradAll(notPolarized),binEdge);

% normalized

bar(binEdge + 0.005, nAligned./n)
axis([binEdge(1)-0.1 binEdge(end) 0 1])
xlabel('dFat')
ylabel('Fraction polarized along dDs');

if saveResults
    saveas(gcf, fullfile(datapath,['alignedFractionHistogramNew_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedFractionHistogramNew_res' num2str(cutoffres) '.fig']));
end
%%
% not normalized

bar(binEdge + 0.005, cat(1,nAligned,nNotAligned,nNotPolarized)', 'stacked')
axis([binEdge(1)-0.1 binEdge(end) 0 max(nAligned)+max(nNotAligned)+2])
xlabel('dFat')
ylabel('number of boundaries');
legend({'Ds aligned','Ds opposed','not polarized'})

if saveResults
    saveas(gcf, fullfile(datapath,['alignedNumberHistogramNew_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedNumberHistogramNew_res' num2str(cutoffres) '.fig']));
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
    saveas(gcf, fullfile(datapath,['alignedFractionHistogramCompatibleNew_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedFractionHistogramCompatibleNew_res' num2str(cutoffres) '.fig']));
end
%%
% not normalized

bar(binEdge + 0.005, cat(1,nAligned,nNotAligned,nNotPolarized)', 'stacked')
axis([binEdge(1) binEdge(end) 0 3])
xlabel('dFat')
ylabel('number of boundaries');
legend({'Ds aligned','Ds opposed','not polarized'})

if saveResults
    saveas(gcf, fullfile(datapath,['alignedNumberHistogramCompatibleNew_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedNumberHistogramCompatibleNew_res' num2str(cutoffres) '.fig']));
end
%%
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
    saveas(gcf, fullfile(datapath,['alignedFractionHistogramIncompatibleNew_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedFractionHistogramIncompatibleNew_res' num2str(cutoffres) '.fig']));
end
%%
% not normalized

bar(binEdge + 0.005, cat(1,nAligned,nNotAligned,nNotPolarized)', 'stacked')
axis([binEdge(1) binEdge(end) 0 20])
xlabel('dFat')
ylabel('number of boundaries');
legend({'Ds aligned','Ds opposed','not polarized'})

if saveResults
    saveas(gcf, fullfile(datapath,['alignedNumberHistogramIncompatibleNew_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedNumberHistogramIncompatibleNew_res' num2str(cutoffres) '.fig']));
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
    saveas(gcf, fullfile(datapath,['alignedFractionCollapsedNew_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedFractionCollapsedNew_res' num2str(cutoffres) '.fig']));
end

% not normalized

bar(cat(1,[nAlignedC nAlignedI],[nNotAlignedC nNotAlignedI],[nNotPolarizedC nNotPolarizedI])', 'stacked');
%axis([0.5 3 0 14])
ylabel('Number of boundaries');
xlabel('Gradients');
set(gca,'XTickLabel',{'compatible','incompatible'})
legend({'Ds aligned','Ds opposed','not polarized'},'Location','NorthWest')

if saveResults
    saveas(gcf, fullfile(datapath,['alignedNumberCollapsedNew_res' num2str(cutoffres) '.png']));
    saveas(gcf, fullfile(datapath,['alignedNumberCollapsedNew_res' num2str(cutoffres) '.fig']));
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
