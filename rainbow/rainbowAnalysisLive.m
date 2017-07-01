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
addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Users/idse/Dropbox/Sprinzak/shared/monoculture movies with Airyscan';
cd(dataDir);

% SETTINGS! 
%-----------------

% % resolution, nm per pixel
% res = 60;

% SETTINGS
thickness = 45;

% emissions peak wavelengths, 530, 610

DsC = 1;
FatC = 2;
DAPIC = [];
files = {...
        fullfile(dataDir, 'SUM_Dual_Dox_post35min(2)_out_Subset1_t15.tif'),...
        fullfile(dataDir, 'SUM_Dual_Dox_post35min(2)_out_Subset1_t20.tif'),...
        fullfile(dataDir, 'SUM_Dual_Dox_post35min(2)_out_Subset1_t25.tif'),...
        fullfile(dataDir, 'SUM_Dual_Dox_post35min(2)_out_Subset1_t30.tif'),...
        fullfile(dataDir, 'SUM_Dual_Dox_post35min(2)_out_Subset1_t35.tif'),...
        fullfile(dataDir, 'SUM_Dual_Dox_post35min(2)_out_Subset1_t40.tif'),...
       };

%% save RGB previews

allRes = {};
Rlim = [];
Glim = [];
Rtol = 0.005;
Gtol = 0.001;
% Rtol = 0.01; above vals are good for rainbows, these vals look better for
% whole image
% Gtol = 0.003;

separations = {};

for fi = numel(files):-1:1 % reverse order to adjust LUT to final

    %-----------------------------
    % read the data
    %-----------------------------

    [filepath,barefname,ext] = fileparts(files{fi});
    [data meta] = readStack(files{fi});
    data = squeeze(data);
    
    dsMIP = data(:,:,DsC);
    fatMIP = data(:,:,FatC);
    
    if fi == 6
        Rlim = stretchlim(mat2gray(dsMIP), Rtol);
        Glim = stretchlim(mat2gray(fatMIP), Gtol);
        Rlim = double(Rlim*(max(dsMIP(:))-min(dsMIP(:))) + min(dsMIP(:)));
        Glim = double(Glim*(max(fatMIP(:))-min(fatMIP(:))) + min(fatMIP(:)));
    end
    
    dsMIP = mat2gray(dsMIP,Rlim);
    fatMIP = mat2gray(fatMIP,Glim);

    RGB = cat(3,dsMIP,fatMIP,0*fatMIP);
    imwrite(RGB, fullfile(dataDir,['RGB_' barefname '.tif']));

    res = round(1000*meta.xres); % resolution in nm/pixel

    %-----------------------------
    % read snake
    %-----------------------------

    %snakeFile = uigetfile('*'); % ['snakes' barefname]
    %snakeFile = fullfile(filepath, snakeFile);

    snakeFile = fullfile(filepath, [barefname '_snakes.txt']);
    if ~exist(snakeFile,'file')
        snakeFile = fullfile(filepath, 'snakes.txt');
    end

    snakes = readSnake(snakeFile, 1);
    nSnakes = size(snakes,2);

    %---------------------------------------------------------
    % take the bow out of rainbow and quantify gradient
    %---------------------------------------------------------
    
    rainbow = {};
    dsCyt = {};
    fatCyt = {};
    nx = {};
    ny = {};

    for i = 1:nSnakes

        % for the rainbow, read out one slice
        d = 60;
        snake = snakes{1, i};
        [fatBdry,~,~,nx{i},ny{i}] = broadImprofile(data(:,:,FatC), snake(:,1), snake(:,2), thickness + 2*d);
        dsBdry = broadImprofile(data(:,:,DsC), snake(:,1), snake(:,2), thickness + 2*d);

        rainbow{i} = cat(3, mat2gray(dsBdry(d+1:end-d,:),Rlim),...
                                mat2gray(fatBdry(d+1:end-d,:),Glim), 0*fatBdry(d+1:end-d,:));

        fth = 300;
        dth = 300;

        r= 8;
        fatMask = imdilate(fatBdry > fth, strel('disk',r));
        dsMask = imdilate(dsBdry > dth, strel('disk',r));
        mask = ~(fatMask | dsMask);
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
    
    clf 
    hold on

    N = 4;
    sep = zeros([1 size(rainbow{1},2)]);

    for y = 1:size(rainbow{1},2)-N

        FatProf = mean(rainbow{i}(:,y:y+N, 1),2);
        DsProf = mean(rainbow{i}(:,y:y+N, 2),2);

        plot(FatProf,'g','LineWidth',3)
        plot(DsProf,'r','LineWidth',1)

        %[maxFatI, maxFatidx] = findpeaks(FatProf./max(FatProf(:)),'MinPeakProminence',0.1,'MinPeakHeight',0.6);
        %[maxDsI, maxDsidx] = findpeaks(DsProf./max(DsProf(:)),'MinPeakProminence',0.1,'MinPeakHeight',0.6);
        [~,maxFatidx] = max(FatProf);
        [~,maxDsidx] = max(DsProf);
        if numel(maxFatidx) == 1 && numel(maxDsidx) == 1
            dmax = -(maxDsidx-maxFatidx);
            [maxDsidx, maxFatidx, dmax];
            sep(y) = dmax;
        else
            sep(y) = NaN;
        end
    end
    hold off

    separations{fi} = sep;
    
    imwrite(rainbow{1}, fullfile(dataDir,['rainbow_' barefname '.tif']));
end

%%
clf
hold on;
N = 1;
d = 50;
ymin = 150;
for y = ymin:ymin+d

    FatProf = mean(rainbow{i}(:,y:y+N, 1),2);
    DsProf = mean(rainbow{i}(:,y:y+N, 2),2);

    xval = (1:numel(FatProf))*res;
    %plot(imopen(FatProf,[1 1 1 1 1]'))
    plot(xval, FatProf,'g')
    plot(xval, DsProf,'r','LineWidth',1)
end
xlabel('distance (nm)');
ylabel('intensity (au)');
set(gcf, 'Color', 'w');
set(gca,'FontSize', 14)
set(gca,'FontWeight', 'bold')
set(gca, 'LineWidth', 2);
hold off;
saveas(gcf,fullfile(dataDir,'xsectionT40_y150_200.png'));
saveas(gcf,fullfile(dataDir,'xsectionT40_y150_200.fig'));

%%

for fi = 1:numel(files)
    
    sep = separations{fi};
    sep(sep > 10) = NaN;
    binedges = -2:10;
    binedgesval = binedges*res;
    n = histc(sep, binedges);
    n = n./sum(n);
    %n = cumsum(n./sum(n));
    [x,y] = histForBarlikePlot(binedges, n');
    eps = 0;
    plot(x*res - fi*eps,y + fi*eps, 'LineWidth',2)
    
    set(gcf, 'Color', 'w');
    set(gca,'FontSize', 14)
    set(gca,'FontWeight', 'bold')
    set(gca, 'LineWidth', 2);
    xlabel('peak separation (nm)')
    ylabel('frequency');

    saveas(gcf,fullfile(dataDir,['separationAll_T' num2str(10+fi*5) '.png']));
    saveas(gcf,fullfile(dataDir,['separationAll_T' num2str(10+fi*5) '.fig']));
end




