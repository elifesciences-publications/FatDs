%--------------------------------------------------------------------------
%
%   script for analyzing FRAP data
%
%   this one is for NCad movies that are like Fat/Ds, so a one sided drop
%   not a window with two sides
%
%--------------------------------------------------------------------------

clear all;
close all;

% add subfolders of folder containing this script to path
[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

cd(testScriptPath);
cd('..');
addpath(genpath(pwd));

%FRAPdataLocations
FRAPdataLocations

autopilot = false;
for movieIdx = 11 

close all;
disp(movies{movieIdx});
if ~exist(fullfile(datapath, movies{movieIdx}),'file')
    disp('does not exist, moving on');
else

[filepath, barefname, ext] = fileparts(fullfile(datapath, movies{movieIdx}));
fname = [barefname ext];
snakeFile = fullfile(filepath, snakes{movieIdx});
cd(datapath);

% % select files 
% [fname, filepath] = uigetfile('*.tif', 'select data file');
%
% cd(filepath)
% snakeFile = uigetfile('*', 'select snake file'); % ['snakes' barefname]
% snakeFile = fullfile(filepath, snakeFile);
% [~,barefname,~] = fileparts(fname);

%-----------------------------
% parameters
%-----------------------------

FatC = 2;                           % channel indices
DsC = 3;   
zSlice = zslices{movieIdx};         % z-slice used for analysis
tmax = cutofftime{movieIdx};        % set temporal cutoff 
framerate = framerates{movieIdx};   % frames/sec
xyresolution = xyres{movieIdx};     % micron / pixel

%-----------------------------
% read the data
%-----------------------------

[data nChannels] = readStack(filepath, fname);
nTimePts = size(data,5);

if isnan(tmax)
    tmax = nTimePts;
end

%%
%-----------------------------
% read snake
%-----------------------------

snake = readSnake(snakeFile, nTimePts);
% backup snake
fullSnake = snake;


% visualize a frame
t = 5;
slice = squeeze(data(:,:,zSlice, DsC, t));
FatSl = squeeze(data(:,:,zSlice, FatC, t));

if ~autopilot
    imshow(cat(3, mat2gray(slice), mat2gray(FatSl), 0*mat2gray(slice)));
    %imshow(slice,[]);
    hold on 
    plot(snake{t}(:,1), snake{t}(:,2), 'c')
    hold off
end

%% prevents sudden jumps in the snake length

%snake = backupsnake;
snake = fullSnake;

if ~isnan(dminMax{movieIdx})
    dminmax = dminMax{movieIdx};
else
    dminmax = 5;
end
dminmax = 2;
for t = 1:tmax-1
    
    DM = distmat(snake{t},snake{t+1});
    % minimal distance of points at t+1 to any point at t
    dmin = min(DM); 
    
    if ~isempty(snake{t+1})
        mask = dmin < dminmax;
        ind = find(mask);
        mask(min(ind):max(ind))=true;
        snake{t+1} = snake{t+1}(mask,:);
    else
        disp(['snake empty at t=' num2str(t+1)]);
    end
end

%% inspect the mask

%figure,
t = 91;
slice = squeeze(data(:,:,zSlice, DsC, t));

%figure(1), 
% imshow(imadjust(slice));
% hold on
% plot(fullSnake{t}(:,1), fullSnake{t}(:,2), 'oc')
% plot(snake{t}(:,1), snake{t}(:,2), '.r')
% hold off

%% plot a profile

t = 2;
thickness = 5;
slice = squeeze(data(:,:,zSlice, FatC, t));
[Iprof, cx, cy] = broadImprofile(slice, fullSnake{t}(:,1), fullSnake{t}(:,2),thickness);
Iprof = sum(Iprof,1);

%figure(2), 
%plot(Iprof)


%%
%-----------------------------------
% read out intensity along snake
%-----------------------------------

% set snake in first bleached frame to snake in prebleach frame
% this is because we cannot register these two frames based on Fat later
% (but maybe Ds based is just as good, anyway)
tmpsnake2 = snake{2};
snake{2} = snake{1};

% get profile and proper distance
Iprofile = {};
propDist = {};

% slice median
sliceMedian = zeros([tmax 1]);

% line thickness
thickness = 9;

for t = 1:tmax

    Fslice = double(squeeze(data(:,:,zSlice, FatC, t)));
    Dslice = double(squeeze(data(:,:,zSlice, DsC, t)));
    
    sliceMedian(t) = median(Fslice(:));
    
%     % possible Gaussian smoothing before readout
%     slice = imfilter(slice, fspecial('gauss', 5, 1));

    % snake points have roughly constant spacing, by taking Npts to be the
    % snake length, we almost have a proper distance parametrization of the
    % intensity profile, making it possible to do registration without
    % first interpolating on a proper distance grid
    %
    % Npts = size(snake{t},1);
    % [cx,cy,c] = improfile(slice, snake{t}(:,1), snake{t}(:,2), Npts,'bicubic');
    
    [FtI, cx, cy] = broadImprofile(Fslice, snake{t}(:,1), snake{t}(:,2), thickness);
    
    % FOR THE MOVIE 2 WHICH HAS REGISTRATION ISSUES
    % FOR TIME POINT 6, check out the point spacing
    % THATS WHERE THE PROBLEM IS
    % scatter(cx, cy)
    
    % maximal intensity projection normal to snake
    [Iprofile{t,1}, idx] = max(FtI,[],1);
    
%     linIdx = sub2ind(size(c), idx, 1:size(c,2));
%     cx = cx(linIdx)';
%     cy = cy(linIdx)';
%   The MIP will not have equally space points which then screws up the
%   registration, so I'll just use the brightest displaced snake for
%   positions
    [~,idx] = max(sum(FtI,2));
    cx = cx(idx,:)';
    cy = cy(idx,:)';

    % calculate the proper distance along the boundary
    dcx = cx - circshift(cx, 1);
    dcx = dcx(2:end);

    dcy = cy - circshift(cy, 1);
    dcy = dcy(2:end);
    
    propDist{t} = [0, cumsum(sqrt(dcx.^2 + dcy.^2))'];
    
    % now also read out Ds at the same positions as Fat
    % but don't take the max intensity but the sum
    ID = broadImprofile(Dslice, snake{t}(:,1), snake{t}(:,2), thickness);
    Iprofile{t,2} = sum(ID,1);
end

%%
% visualize a frame
t = 8;

slice = squeeze(data(:,:,zSlice, FatC, t));
Ds = mat2gray(squeeze(data(:,:,zSlice, DsC, t)));

h = figure;
if autopilot 
	set(h,'Visible','off');
end
imshow(cat(3, mat2gray(Ds), mat2gray(slice), 0*Ds),[]);
%figure(1), 
%imshow(slice,[]);
hold on
plot(snake{t}(:,1), snake{t}(:,2), 'c')
hold off
saveas(gcf,fullfile(filepath, ['snakeOnImage_t' num2str(t) '.png']));

%%
%c = medfilt2(c, [3 1]);
%figure(2), 
% figure,
% t = 2;
% plot(propDist{t}, 4*Iprofile{t,1}, 'g')
% hold on;
% plot(propDist{t}, Iprofile{t,2}, 'r')
% hold off;
% legend({'Fat','Ds'})

%%
%-----------------------------------
% register frames
%-----------------------------------

relShift = zeros([1 tmax]);
totalShift = zeros([1 tmax]);
propDistReg = {};
propDistReg{1} = propDist{1};

for t = 1:tmax-1
    
    if t > 1
        
        f1 = Iprofile{t};
        f2 = Iprofile{t+1};
        
        % below block can be removed to register based on intensity
        % like this it registers the derivative over the bulk of the image
        % to prevent the edge of the snake messing up registration
        % on some movies it works better with this, on some better without
%         %------- 
%         sig = 1;
%         df2 = conv(f2,GaussD(sig,1,2),'same');
%         df1 = conv(f1,GaussD(sig,1,2),'same');
%         x = round(3*sig);
%         f2 = df2(x:end-x);
%        	f1 = df1(x:end-x);
%         %-------

        L1 = length(f1);
        L2 = length(f2);
        
        if L1 > L2 
            pad = zeros([1 L1-L2]);
            f2 = [f2 pad];
        elseif L1 < L2
            pad = zeros([1 L2-L1]);
            f1 = [f1 pad];
        end
        
        relShift(t) = register1D(f1', f2', 101);
        totalShift(t) = totalShift(t-1) + relShift(t);
    end
    
    % in principle ds (proper distance differential, not dachsous)
    % is not constant so this is a little sloppy, but a very
    %good approximation
    ds = propDist{t+1}(end)/length(propDist{t+1});
    propDistReg{t+1} = propDist{t+1} + totalShift(t)*ds;
end

% %%
% % ------- check result ------- 
% t = 6;
% dt = 1;
% 
% c = 2;
% f1 = Iprofile{t,c};
% f2 = Iprofile{t+dt,c};
% 
% %plot(1:length(f1), f1, 1:length(f2), f2)
% 
% plot(propDistReg{t}, f1, 'b')
% hold on
% plot(propDistReg{t+dt},f2, 'r')
% plot(propDist{t+dt}, f2, '--r')
% hold off
% %axis([0 50  3000 20000])

%%
%---------------------------------------------------
% spacetime profile
%---------------------------------------------------

% an interpolation step is required to get the registered frames on a
% square grid

smin = 0;
smax = 0;
for t = 1:tmax
   smin = min(smin, min(propDistReg{t}));
   smax = max(smax, max(propDistReg{t}));
end
smin = floor(smin);
smax = ceil(smax);
% ds = (smax-smin)/(2*Npts);
ds = 1;

[S,T] = meshgrid(smin:ds:smax, 1:tmax);
FtI = zeros(size(S));
FtIraw = zeros(size(S));
DsIraw = zeros(size(S));

medfiltwidth = 2;

for t = 1:tmax
    
    FtIraw(t,:) = interp1(propDistReg{t}, Iprofile{t,1}, S(t,:), 'linear');
    FtIraw(isnan(FtIraw)) = 0;
    
    % Dachsous too
    DsIraw(t,:) = interp1(propDistReg{t}, Iprofile{t,2}, S(t,:), 'linear');
    DsIraw(isnan(DsIraw)) = 0;
end

% %----------------------------------------
% % 3D plot to visualize temporal evolution
% %----------------------------------------
% 
% surf(S,T,FtIraw)
% view([1 1 1])
% view([0 0 -1])
% axis([smin smax 1 tmax 1 max(FtIraw(:))]);
% colormap jet

imshow(FtIraw,[])

% save for later inspection
imwrite(mat2gray(FtIraw), fullfile(filepath,['profile_' barefname '_raw.tif']));

% %%
% %-----------------------------
% % read the FRAP location
% %-----------------------------
% 
% w = 5;
% 
% rgnfname = [barefname '.rgn'];
% [xF yF wF hF] = readFRAPlocation(filepath, rgnfname);
% 
% % part of the curve inside the FRAP region
% t=1;
% FRAPpoly = [[xF xF xF+w xF+wF]', [yF yF+hF yF+hF yF]'];
% inFRAP = inpolygon(snake{t}(:,1), snake{t}(:,2), FRAPpoly(:,1), FRAPpoly(:,2));
% 
% % the proper distance of the FRAP edge
% sFRAP = [min(propDist{t}(inFRAP)) max(propDist{t}(inFRAP))];
% 
% % the edge of the FRAP region that is closest to the middle is the one we
% % want
% [~,mi] = min(abs((sFRAP - smin) - (smax - sFRAP)));
% 
% %figure
% % the FRAP edge doesn't lie right on the brightness edge
% 
% % --- PAY ATTENTION HERE ----
% % FI , FIraw, FIcorr
% % imshow(FtIraw,[])
% % hold on;
% % plot((sFRAP(mi)-smin+1)*ones([1 tmax])/ds, 1:tmax, 'y');
% % hold off;

%%
%-----------------------------
% fit erf to space
%-----------------------------

% % stupid stuff for Fat-Ds -> move to parameter file
% if movieIdx == 6
%     idx = 68; 
%     sgn = 1;
% elseif movieIdx == 8
%     idx = 158;
%     sgn = -1;
% elseif movieIdx == 9
%     idx = 58;
%     sgn = -1;
% elseif movieIdx == 11mo
%     idx = 41;
%     sgn = -1;
% else
%     sgn = 1;
%     [~,cidx] = find(S > sFRAP(mi));
%     idx = min(cidx);
% end

if ~isnan(sign{movieIdx})
    sgn = sign{movieIdx};
else
    sgn = 1;
end
if ~isnan(windowcenter{movieIdx})
    idx = windowcenter{movieIdx};
else
    idx = floor(size(FtIraw,2)/2);
end
if ~isnan(windowwidth{movieIdx})
    L = windowwidth{movieIdx};
else
    L = 40;
end

center = idx;%round(mean(find(FtIraw(tidx,:) > 0)));
lowerIdx = max(1,center - L);
upperIdx = min(center + L, size(FtIraw,2));

xidx = lowerIdx:upperIdx;
x = xidx*xyresolution;

% function f to be fitted
%
% f = U0*(1 - A*[1 + erf((x-x0)/L)]/2 )
%
% where L = sqrt(4 D t) for a sharp initial condition
% and L = sqrt(1 + 4 m^2 D t)/m for a profile of width 1/m
% A = b exp(-gamma t)
%
% p = [U0 A x0 L]
f = @(p,x) p(1)*(1 - p(2)*(1 + erf((x - p(3))/p(4)))/2);

param = struct('U0', [], 'A', [], 'x0', [], 'L', [], 'D', []);

tmaxFit = cutofftime{movieIdx};
badones = false([1 tmaxFit-1]);
exitflag = 1;

for tidx = 2:tmaxFit

    % a hack to cutoff the data beyond the plateau
    if sgn < 0
        minidx = 1;
        maxidx = xidx(1) - 1 + find(imextendedmax(FtIraw(tidx,xidx),...
                            max(FtIraw(tidx,xidx))/2),1,'last');
    else
        minidx = xidx(1) - 1 + find(imextendedmax(FtIraw(tidx,xidx),...
                            max(FtIraw(tidx,xidx))/2),1,'first');
        maxidx = max(xidx);
    end
    
    % cutoff x where there is no signal
    goodx = FtIraw(tidx,xidx)>0 & xidx < maxidx & xidx > minidx; 
    xidxt = xidx(goodx);
    xt = x(goodx);

    % E: energy functional = sum(fun.^2)
    fun = @(p) f(p,xt) - FtIraw(tidx,xidxt);

    if tidx==2 || exitflag ==0
        U00 = max(FtIraw(tidx,xidxt));
        A0 = 0.8; % 1
        x0 = mean(xt);
        L0 = sgn*2; % 5
    else 
        U00 = param(tidx-1).U0;
        A0 = param(tidx-1).A;
        x0 = param(tidx-1).x0;
        L0 = param(tidx-1).L;
    end
  
    pinit = [U00 A0 x0 L0];

    % lower and upper bounds
    lb = -Inf + pinit*0;
    ub = Inf + pinit*0;
    lb(2) = 0;
    ub(2) = 1;
    ub(1) = 2*max(FtIraw(:));
    lb(1) = min(FtIraw(:));
    
    %-----
    options.Algorithm = 'trust-region-reflective';
    options.TolFun = 1e-1;
    options.Display = 'off';
    
    if numel(xidxt) > 10
        p = lsqnonlin(fun, pinit, lb, ub,options);
    else
        p = NaN*pinit;
        exitflag = 0;
    end

    if exitflag == 0
        badones(tidx-1) = true;
    end
    
    param(tidx).U0  = p(1);
    param(tidx).A  = p(2);
    param(tidx).x0  = p(3);
    param(tidx).L = p(4);
    
    param(tidx).D = p(4)^2/(4*(tidx/framerate));
end

%%
for tidx=2:tmaxFit
    h = figure;
    if autopilot 
        set(h,'Visible','off');
    end
    plot(x, FtIraw(tidx,xidx),'--b')
    %axis([x(1) x(end) 500 3000]);
    hold on
    p = [param(tidx).U0, param(tidx).A, param(tidx).x0, param(tidx).L];
    plot(x, f(p,x), 'r');
    %plot(x,(xidx < 50)*10^4)
    hold off
    saveas(h, fullfile(filepath,['fit_t' num2str(tidx) '.png']));
    close;
end

%%
% FRAP boundary width 
t = (2:tmaxFit)/framerate;
plot(t,sgn*[param(:).L])
axis([0 t(end) 0 5]);
xlabel('time (seconds)')
ylabel('boundary width (micron)');
%saveas(gcf, fullfile(filepath, 'bdryWidthVsTime.png'));

%%
% FRAP amplitude
t = (2:tmaxFit)/framerate;
plot(t,[param(:).A])
%axis([0 t(end) 0 20]);
xlabel('time (seconds)')
ylabel('amplitude');
saveas(gcf, fullfile(filepath, 'amplitudeVsTime.png'));

%%
% log FRAP amplitude
t = (2:tmaxFit)/framerate;
plot(t,log([param(:).A]))
%axis([0 t(end) 0 20]);
xlabel('time (seconds)')
ylabel('log(amplitude)');
saveas(gcf, fullfile(filepath, 'logAmplitudeVsTime.png'));

%% fit amplitude for off rate

t = (2:tmaxFit)/framerate;
logA = log([param.A]);
logA = logA(1:tmaxFit-1);

logA = logA(~badones);
t = t(~badones);

N = length([param.L]);

alpha = sum(t.^2);
beta = N;
gamma = sum(t);

M = [alpha, gamma; gamma, beta];

Minv = inv(M);

p = sum(t.*logA);
q = sum(logA);

v = [p q]';

paramFit = Minv*v;

gamma = -paramFit(1);
logb = paramFit(2);

% error CHECK DATA ANALYSIS BOOK FOR FACTORS 2
sig = sqrt(sum((logA - logb + gamma*t).^2)/(N-1));
cov = 2*Minv*sig^2;
sigGamma = sqrt(cov(1,1));
sigLogb = sqrt(cov(2,2));

allParam = struct();
allParam.gamma = [gamma sigGamma];
allParam.logb = [logb sigLogb];

% visualize
h = figure;
if autopilot 
	set(h,'Visible','off');
end
plot(t,logA)
hold on
plot(t,logb - gamma*t,'--r')
hold off
axis([1 t(end) min(logA(:))-0.05 max(logA(:))+0.05])
xlabel('time (seconds)')
ylabel('amplitude');
legend({'logA', 'logA = logb - \gamma t'});
title(['log b = ' num2str(logb,'%.2e') ', \gamma = ' num2str(gamma,'%.2e') '(' num2str(sigGamma,'%.2e') ')']);
saveas(gcf, fullfile(filepath, 'amplitudeVsTime.png'));

%%
% the profile steepness doesn't change over time

% and L = sqrt(1 + 4 m^2 D t)/m for a profile of initial width 1/m = L0
% L^2 = L0^2 + 4 D t
% so we're fitting a straight line

% see book on LSQ, there m = 4D, c = L0^2, Y = L^2

t = (2:tmaxFit)/framerate;
Lsq = [param.L].^2;
Lsq = Lsq(1:tmaxFit-1);

Lsq = Lsq(~badones);
t = t(~badones);

N = length([param.L]);

alpha = sum(t.^2);
beta = N;
gamma = sum(t);

M = [alpha, gamma; gamma, beta];

Minv = inv(M);

p = sum(t.*Lsq);
q = sum(Lsq);

v = [p q]';

paramFit = Minv*v;

Draw = paramFit(1)/4;
L0 = sqrt(paramFit(2));

% take into account the prior
D = max(0, Draw);
D
sig = sqrt(sum((Lsq - L0^2 - 4*D*t).^2)/(N-1));
sigraw = sqrt(sum((Lsq - L0^2 - 4*Draw*t).^2)/(N-1));

cov = 2*Minv*sig^2;
covraw = 2*Minv*sigraw^2;

sigDraw = sqrt(covraw(1,1))
sigD = sqrt(cov(1,1))
sigL = sqrt(cov(2,2));

allParam.D = [D sigD];
allParam.Draw = [Draw sigDraw];
allParam.L0 = [L0 sigL];

% visualize
h = figure;
if autopilot 
	set(h,'Visible','off');
end
plot(t,Lsq)
hold on
plot(t,L0^2 + 4*Draw*t,'--r')
hold off
axis([1 t(end) 0 max(Lsq(:))])
xlabel('time (seconds)')
ylabel('square width (micron^2)');
legend({'L^2', 'L^2 = L0^2 + 4 D t'});
title(['L0^2 = ' num2str(L0^2,'%.2e') ', D = ' num2str(Draw,'%.2e') '(' num2str(sigDraw,'%.2e') ')']);
saveas(gcf, fullfile(filepath, 'squareWidthVsTime.png'));

% save fitted parameters
save(fullfile(filepath,'result.mat'),'param','allParam','FtIraw');

end
end