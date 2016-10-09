%--------------------------------------------------------------------------
%
%   script for analyzing FRAP data
%
%--------------------------------------------------------------------------

clear all;
close all;

% add subfolders of folder containing this script to path
[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

cd(testScriptPath);
cd('..');
addpath(genpath(pwd));

datapath = '/Users/idse/Dropbox/Sprinzak/data/';
%datapath = '/Users/idse/Dropbox/Sprinzak/shared/frap movies with snakes/2/';
cd(datapath);

% select files 
[fname, filepath] = uigetfile('*.tif', 'select data file');

cd(filepath)
snakeFile = uigetfile('*', 'select snake file'); % ['snakes' barefname]
snakeFile = fullfile(filepath, snakeFile);

[~,barefname,~] = fileparts(fname);

% channel indices
FatC = 2;
DsC = 3;
%FatC = 2;
%DsC = 1;

%-----------------------------
% read the data
%-----------------------------

[data nChannels] = readStack(filepath, fname);
nTimePts = size(data,5);

%%
%-----------------------------
% OLGA! SETTINGS
%-----------------------------

% set temporal cutoff 
tmax = 80;%nTimePts;
framerate = 1/5; % frames/sec
xyresolution = 0.27; % micron / pixel
zSlice = 1;

%-----------------------------
% read snake
%-----------------------------

snake = readSnake(snakeFile, nTimePts);

%%
% visualize a frame
t = 7;
slice = squeeze(data(:,:,zSlice, DsC, t));
FatSl = squeeze(data(:,:,zSlice, FatC, t));

imshow(cat(3, mat2gray(slice), mat2gray(FatSl), 0*mat2gray(slice)));
%imshow(slice,[]);
hold on 
plot(snake{t}(:,1), snake{t}(:,2), 'c')
hold off

% backup snake
fullSnake = snake;

%%
%-----------------------------------
% crop snake based on unbleached Ds
%-----------------------------------

nsig = 2;
openSize = 6;
dilateSize = openSize;

for t = 1:tmax
    
    snake{t} = fullSnake{t};
    slice = double(squeeze(data(:,:,zSlice, DsC, t)));

    thresh = mean(slice(:)) + nsig*std(slice(:));

    Iprof = broadImprofile(slice, snake{t}(:,1), snake{t}(:,2));
    snakeMask = Iprof > thresh;
    snakeMask = imopen(snakeMask', true([1 openSize])');
    snakeMask = imdilate(snakeMask, true([1 dilateSize])');

    newSnake = snake{t};
    
    if all(snakeMask)

        newSnake = growSnake(newSnake, 5);
        
        Iprof = broadImprofile(slice, newSnake(:,1), newSnake(:,2));
        snakeMask = Iprof > thresh;
        snakeMask = imopen(snakeMask', true([1 openSize])');
        snakeMask = imdilate(snakeMask, true([1 dilateSize])');  
    end
    
    snake{t} = newSnake(snakeMask,:);
end

% set snake in first bleached frame to snake in prebleach frame
% this is because we cannot register these two frames based on Fat later
% (but maybe Ds based is just as good, anyway)
snake{2} = snake{1};


%% inspect the mask

t = 6;
slice = squeeze(data(:,:,zSlice, DsC, t));

% %figure(1), 
% imshow(slice,[]);
% hold on
% plot(fullSnake{t}(:,1), fullSnake{t}(:,2), 'oc')
% plot(snake{t}(:,1), snake{t}(:,2), '.r')
% hold off

[Iprof, cx, cy] = broadImprofile(slice, fullSnake{t}(:,1), fullSnake{t}(:,2));
snakeMask = Iprof > thresh;
snakeMask = imopen(snakeMask', true([1 openSize])');
snakeMask = imdilate(snakeMask, true([1 dilateSize])');

%figure(2), 
plot(Iprof)
hold on;
plot(1000*snakeMask, 'g');
hold off;

% %%
% %-----------------------------------
% % try out my function broadImprofile
% %-----------------------------------
% 
% % fat as in thickness, not the protein
% 
% thickness = 7;
% [c,cx,cy] = broadImprofile(slice, snake{t}(:,1), snake{t}(:,2), thickness);
% 
% figure(1), imshow(slice,[]);
% hold on
% plot(cx', cy', 'r')
% hold off
% 
% figure(2)
% imshow(c,[])

%%
%-----------------------------
% make a film strip
%-----------------------------

% visualize a frame
stripTimes = [1 2 40]; 
xmin = 165;
xmax = xmin + 80;
ymin = 325;
ymax = ymin + 80;
gammaval = 0.5;
s = 10;

for t= stripTimes
    
    h = figure;%('Position',[1 1 s*(xmax-xmin+1) s*(ymax-ymin+1)])
    
    slice = imadjust(squeeze(data(ymin:ymax,xmin:xmax,zSlice, DsC, t)), [0 1], [0 1],gammaval);
    FatSl = imadjust(squeeze(data(ymin:ymax,xmin:xmax,zSlice, FatC, t)), [0 1], [0 1],gammaval);
    imshow(cat(3, mat2gray(slice), mat2gray(FatSl), 0*mat2gray(slice)),...
                            'InitialMagnification', 300);

    hold on 

    % box around accumulating interface
    thickness = 9;
    [~,cx,cy] = broadImprofile(slice, snake{t}(:,1)- xmin, snake{t}(:,2)- ymin, thickness);
    plot([cx(1,:) cx(end,end:-1:1) cx(1,1)], [cy(1,:) cy(end,end:-1:1) cy(1,1)], '--c', 'LineWidth',2)

    % time label
    textmargin = 7;
    fsize = 25;
    text(textmargin,  textmargin, ['t=' num2str((t-1)/framerate) 's'], 'Color', 'white','FontSize',fsize,'FontWeight','Bold')

    if t == 1
        % Fat Ds labels
        text(textmargin, ymax - ymin - textmargin, 'Fat4', 'Color', 'green','FontSize',fsize,'FontWeight','Bold')
        text(ymax - ymin - textmargin-5, textmargin, 'Ds1', 'Color', 'red','FontSize',fsize,'FontWeight','Bold')

        % scalebar
        fivemuinpixels = 5/xyresolution;
        barxval = xmax-xmin-textmargin-fivemuinpixels:0.1:xmax-xmin-textmargin;
        baryval = (ymax-ymin-textmargin)*ones(size(barxval));
        plot(barxval, baryval,'Color','white','LineWidth',2)
        text(barxval(1)+3, baryval(1)+3, '5\mum', 'Color', 'white','FontSize',fsize,'FontWeight','Bold')
    end
    set(gcf,'color','w');
    hold off
    
    saveas(h, [barefname '_film_T' num2str(t) '.png']);
end

% backup snake
fullSnake = snake;

%%
%-----------------------------------
% read out Fat and Ds along snake
%-----------------------------------

% get profile and proper distance
Iprofile = {};
propDist = {};

% slice median
sliceMedian = zeros([tmax 1]);

% line thickness
thickness = 5;

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
    
    [FI, cx, cy] = broadImprofile(Fslice, snake{t}(:,1), snake{t}(:,2), thickness);
    
    % FOR THE MOVIE 2 WHICH HAS REGISTRATION ISSUES
    % FOR TIME POINT 6, check out the point spacing
    % THATS WHERE THE PROBLEM IS
    % scatter(cx, cy)
    
    % maximal intensity projection normal to snake
    [Iprofile{t,1}, idx] = max(FI,[],1);
    
%     linIdx = sub2ind(size(c), idx, 1:size(c,2));
%     cx = cx(linIdx)';
%     cy = cy(linIdx)';
%   The MIP will not have equally space points which then screws up the
%   registration, so I'll just use the brightest displaced snake for
%   positions
    [~,idx] = max(sum(FI,2));
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
t = 6;

slice = squeeze(data(:,:,zSlice, FatC, t));
Ds = mat2gray(squeeze(data(:,:,zSlice, DsC, t)));

%figure(1), 
imshow(cat(3, mat2gray(Ds), mat2gray(slice), 0*Ds),[]);
%figure(1), 
%imshow(slice,[]);
hold on
plot(snake{t}(:,1), snake{t}(:,2), 'c')
hold off

%%
%c = medfilt2(c, [3 1]);
%figure(2), 
t = 4;
plot(propDist{t}, 4*Iprofile{t,1}, 'g')
hold on;
plot(propDist{t}, Iprofile{t,2}, 'r')
plot(propDist{t}, thresh*ones([1 length(propDist{t})]), 'b')
hold off;

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

%%
% ------- check result ------- 
t = 5;
dt = 1;

c = 2;
f1 = Iprofile{t,c};
f2 = Iprofile{t+dt,c};

%plot(1:length(f1), f1, 1:length(f2), f2)

plot(propDistReg{t}, f1, 'b')
hold on
plot(propDistReg{t+dt}, f2, 'r')
plot(propDist{t+dt}, f2, '--r')
hold off
%axis([0 50  3000 20000])

%%
%-----------------------------------
% overall bleach correction
%-----------------------------------

bg = 1.5*sliceMedian;
% there's barely any variation so
bg = 0*bg + mean(bg);

% total image intensity to estimate photobleaching
Itot = zeros([tmax 1]);

for t = 1:tmax
    
    slice = squeeze(data(:,:,zSlice, FatC, t));
    foregroundMask = bwareaopen(slice > bg(t), 4);
    Itot(t) = sum(slice(foregroundMask) - bg(t));
end

%max(Itot(:))/min(Itot(:))
bleachFactor = Itot/max(Itot(:));

% figure, plot(bleachFactor)
% FI = getframe(gcf);
% imwrite(FI.cdata, ['bleachFactor_' barefname '.png']);
% %close

%%
% imshow(mat2gray(slice).*mat2gray(foregroundMask))

% %% manual bg selection to compare
% 
% Fat = squeeze(data(:,:,zSlice, FatC, t));
% Ds = squeeze(data(:,:,zSlice, DsC, t));
% colorSlice = cat(3, mat2gray(Ds), mat2gray(Fat), mat2gray(0*Fat));
% 
% imshow(colorSlice);
% rect = getrect;
% rect = round(rect);
% slicebg = Fat(rect(1):rect(1)+rect(3),rect(2):rect(2)+rect(4));
% mean(slicebg(:))

%%
%---------------------------------------------------
% bleach and background corrected spacetime profile
%---------------------------------------------------
figure
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
FI = zeros(size(S));
FIraw = FI;
DIraw = FI;

medfiltwidth = 2;

for t = 3:tmax
    
    bcorr = 1./bleachFactor(t);
    FIraw(t,:) = interp1(propDistReg{t}, Iprofile{t,1}, S(t,:), 'linear');
    FIraw(isnan(FIraw)) = 0;
    
    % Dachsous too
    DIraw(t,:) = interp1(propDistReg{t}, Iprofile{t,2}, S(t,:), 'linear');
    DIraw(isnan(DIraw)) = 0;
    
    % subtract background and multiply by bleachcorrection
    FI(t,:) = bcorr*max(FIraw(t,:) - bg(t),0);
    
    % then median filter
    FI(t,:) = medfilt2(FI(t,:),[1 medfiltwidth]);
end

%----------------------------------------
% 3D plot to visualize temporal evolution
%----------------------------------------

surf(S,T,FI)
view([1 1 1])
view([0 0 -1])
axis([smin smax 1 tmax 1 max(FI(:))]);
colormap jet

% save for later inspection

imwrite(mat2gray(FIraw), ['profile_' barefname '_raw.tif'])
imwrite(mat2gray(FI), ['profile_' barefname '_bleachcorr.tif'])

%% normalize by other channel to correct for z-drift and bleaching

w = 5;
FIcorr = FIraw;
DIcorr = DIraw;
FIbgsub = FIraw;

Fbg = 250;
Dbg = 300;

for t =1:size(FIraw,1)
    
    FIcorr(t,:) = medfilt2(FIraw(t,:), [1 1]);
    %FIcorr(t,:) = FIraw(t,:);    
    %Fbg = min(FIcorr(t,FIcorr(t,:)>0));
    
    FIcorr(t,:) = FIcorr(t,:) - Fbg;
    FIbgsub(t,:) = FIcorr(t,:);
    
    DIcorr(t,:) = medfilt2(DIraw(t,:), [1 w]);
    %Dbg = min(DIcorr(t,DIcorr(t,:)>0));
    
    DIcorr(t,:) = DIcorr(t,:) - Dbg;
    DIcorr(t,DIcorr(t,:)<0)=0;
    
    FIcorr(t,:) = FIcorr(t,:)./DIcorr(t,:);
    %FIcorr(t,:) = (max(FIraw(t,:)-Fbg)/max(FIcorr(t,:)))*FIcorr(t,:) + Fbg;
end

FIcorr = (max(FIraw(:))-Fbg)/max(FIcorr(:))*FIcorr + Fbg;

%FIcorr = medfilt2(FIcorr, [2 1]);
t = 6;
plot(FIcorr(t,:),'b')
hold on
plot(FIraw(t,:),'--g');
plot(FIbgsub(t,:),'g');
plot(DIraw(t,:),'--r')
plot(DIcorr(t,:),'r')
hold off
axis([0 82 0 8000]);

%%
%-----------------------------
% read the FRAP location
%-----------------------------

rgnfname = [barefname '.rgn'];
[xF yF wF hF] = readFRAPlocation(filepath, rgnfname);

% part of the curve inside the FRAP region
t=1;
FRAPpoly = [[xF xF xF+w xF+wF]', [yF yF+hF yF+hF yF]'];
inFRAP = inpolygon(snake{t}(:,1), snake{t}(:,2), FRAPpoly(:,1), FRAPpoly(:,2));

% the proper distance of the FRAP edge
sFRAP = [min(propDist{t}(inFRAP)) max(propDist{t}(inFRAP))];

% the edge of the FRAP region that is closest to the middle is the one we
% want
[~,mi] = min(abs((sFRAP - smin) - (smax - sFRAP)));

figure
% the FRAP edge doesn't lie right on the brightness edge

% --- PAY ATTENTION HERE ----
% FI , FIraw, FIcorr
imshow(FIraw,[])
hold on;
plot((sFRAP(mi)-smin+1)*ones([1 tmax])/ds, 1:tmax, 'y');
hold off;

%%
%-----------------------------
% kymograph figure
%-----------------------------

% FIcorrfilt = medfilt2(FIcorr,[2 1],'symmetric');
% FIcorrfilt(1,:) = FIcorr(1,:);
h = figure;
tmaxKymo = 70;
FIcorrfilt = FIcorr(1:tmaxKymo,:);

scale = 10;
%viewvec = [1 1 1];
viewvec = [0 0 1];

xidx = 1:size(FIcorr,2);
x = xidx*xyresolution;
tbleach = 1;
tidx = 1:tmaxKymo;
t = (tidx - tbleach)/framerate;

[X,T] = meshgrid(x-min(x),t);
Z = FIcorrfilt;
surf(X,T,Z, 'EdgeColor','none')
axis equal;
set(gca, 'DataAspectRatio', [1 16 1]);
xlabel('position (\mum)','FontSize',fsize);
ylabel('time (s)','FontSize',fsize);
view(viewvec)
%title('Kymograph','FontSize',fsize);
% to get the same colorscale
Zmin = min(Z(:));
Zmax = max(Z(:));
colormap(gray(256));
set(gcf,'color','w');
set(gca,'FontSize', fsize)
set(gca,'FontWeight', 'bold')
set(gca, 'LineWidth', 2);
grid off;

saveas(h, [barefname '_kymograph.tif']);

%%
% visualize in slice
tidx = 2;
slice = squeeze(data(:,:,zSlice, FatC, tidx));
%figure(1), 
imshow(slice,[])
rectangle('Position',[xF, yF, wF, hF], 'EdgeColor', 'y')
hold on
plot(snake{tidx}(:,1), snake{tidx}(:,2), 'g');
hold off

% plot frapped part of curve in red
hold on
plot(snake{tidx}(inFRAP,1), snake{tidx}(inFRAP,2), 'r');
hold off

% get the profile
%slice = imfilter(slice, fspecial('gauss', 10, 2)); 
Npts = length(snake{tidx});
FIframe = improfile(slice, snake{tidx}(:,1), snake{tidx}(:,2), Npts, 'bicubic');
cFRAP = FIframe;
cFRAP(~inFRAP) = NaN;

% % plot intensity profile
% figure(2), plot(c, 'g')
% hold on
% plot(cFRAP, 'r');
% hold off

%%
%-----------------------------
% overall recovery
%-----------------------------

% get total intensity in part of the FRAP window
w = 12;

% margin away from the FRAP edge
marg = 5;

% idx of the FRAP edge in S
[~,cidx] = find(S > sFRAP(mi));
idx = min(cidx);
nmi = setdiff([1 2], mi);
FRAPsign = sign(sFRAP(nmi) - sFRAP(mi));

FRAPmask = false(size(S));
FRAPrange = idx + FRAPsign*marg:FRAPsign:idx + FRAPsign*(w+marg);
FRAPmask(:, FRAPrange) = true;

notFRAPmask = false(size(S));
notFRAPrange = idx - FRAPsign*marg:-FRAPsign:idx - FRAPsign*(w + marg);
notFRAPmask(:, notFRAPrange) = true;

%figure(1), 
imshow(FI.*(FRAPmask | notFRAPmask), [])
hold on;
plot((sFRAP(mi)-smin+1)*ones([1 tmax])/ds, 1:tmax, 'y');
hold off;

% subtract background
Ifrapped = FI.*FRAPmask;
Ifrapped(isnan(Ifrapped))=0;
IfrappedTot = sum(Ifrapped, 2);

InotFrapped = FI.*notFRAPmask;
InotFrapped(isnan(InotFrapped)) = 0;
InotFrappedTot = sum(InotFrapped,2);

% figure(1)
% imshow(I,[])
% colors = lines(3);
% hold on;
% plot((sFRAP-smin+1)*ones([1 tmax])/ds, 1:tmax, 'Color', colors(1,:));
% plot((sFRAP + w -smin+1)*ones([1 tmax])/ds, 1:tmax, 'Color', colors(1,:));
% plot(idx - FRAPsign*marg*ones([1 tmax]), 1:tmax, 'Color', colors(2,:));
% plot((sFRAP - w -smin+1)*ones([1 tmax])/ds, 1:tmax, 'Color', colors(2,:));
% hold off;

%%
% fitting
%----------------

tbleach = 1;
tdata = (1:tmax)/framerate;
Idata = IfrappedTot(1:tmax)/IfrappedTot(1);

% function f to be fitted
% parameters p = [A tau]
f = @(p,t) p(1)*(1-exp(- (t'-tbleach)/p(2) ));

% E: energy functional
E = @(p) sum((Idata(tbleach+1:end) - f(p,tdata(tbleach+1:end))).^2);

% initial estimates for parameters
pinit = [0.5, 300];

%optimise
[p,fminres] = fminsearch(E, pinit);

% plot
figure(2)

IdataNotFrapped = InotFrappedTot(1:tmax)/IfrappedTot(1);
%IdataNotFrapped = IdataNotFrapped./bleachFactor;

plot(tdata, Idata, tdata, IdataNotFrapped, tdata, f(p,tdata));

% hold on;
% plot(tdata,bleachFactor);
% hold off;

xlabel('time (seconds)');
ylabel('intensity');

legend( 'Bdry Fat in bleach ROI',...
        'Bdry Fat not bleached',...
        ['A(1-e^{-t/\tau}), A=' num2str(p(1),2) ', \tau=' num2str(p(2),3)]);

% SUBSTRATE DEPENDENT RECOVERY?

%%
%-----------------------------
% fit erf to space
%-----------------------------

lowerIdx = min([FRAPrange notFRAPrange]);
upperIdx = max([FRAPrange notFRAPrange]);


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

for tidx = 2:tmax

    % E: energy functional
    E = @(p) sum((FI(tidx,xidx) - f(p,x)).^2);

    U00 = max(FI(tidx,:));
    A0 = 1;
    x0 = mean(x);
    m0 = 5;

    pinit = [U00 A0 x0 m0];

    [p,fminres,exitflag] = fminsearch(E, pinit);

    if exitflag == 0
        tidx
    end
    
    param(tidx).U0  = p(1);
    param(tidx).A  = p(2);
    param(tidx).x0  = p(3);
    param(tidx).L = p(4);
    
    param(tidx).D = p(4)^2/(4*(tidx/framerate));
end

%%
tidx=70;
plot(x, FI(tidx,xidx))
hold on
plot(x, FIraw(tidx,xidx),'--b')
plot(x, FIcorr(tidx,xidx),'--g')
p = [param(tidx).U0, param(tidx).A, param(tidx).x0, param(tidx).L];
plot(x, f(p,x), 'r')
hold off

%%
% the profile steepness doesn't change over time

% see book, m = 4D, c = L0^2, Y = L^2

tmaxFit = 30;
t = (2:tmaxFit)/framerate;
Lsq = [param.L].^2;
Lsq = Lsq(1:(tmaxFit-1));

plot(Lsq)
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

D = paramFit(1)/4
L0 = sqrt(paramFit(2))

sig = sqrt(sum((Lsq - L0^2 - 4*D*t).^2)/(N-1));

cov = 2*Minv*sig^2;
sigD = sqrt(cov(1,1));
sigL = sqrt(cov(2,2));

allParam = struct();
allParam.DperTime = [D sigD];
allParam.L0perTime = [L0 sigL];

%%
%----------------------------------------
% fit recovery profile in space and time
%----------------------------------------

lowerIdx = min([FRAPrange notFRAPrange]);
upperIdx = max([FRAPrange notFRAPrange]);

xidx = lowerIdx:upperIdx;
x = xidx*xyresolution;
tbleach = 1;
tidx = 2:tmax;
t = (tidx - tbleach)/framerate;

% function f to be fitted
%
% f = U0*(1 - b*exp(-gamma t)*[1 + erf((x-x0)/L)]/2 )
%
% L is a constant here, whereas it should be L = sqrt(L0^2 + 4 D t)
% which is the next step
%
% p = [U0 b gamma x0 L]
expfactor = @(p,x,t) p(2)*exp(-p(3)*t)'/2;
f = @(p,x,t) p(1)*(1 - expfactor(p,x,t)*(1 + erf((x - p(4))/p(5))));

% E: energy functional
E = @(p) sum(sum((FIcorr(tidx,xidx) - f(p,x,t)).^2));

U00 = max(FI(:));
b0 = 1;
gamma0 = 1/300;
x0 = mean(x);
m0 = 1;

pinit = [U00 b0 gamma0 x0 m0];

[p,fminres] = fminsearch(E, pinit);

% visualize result

scale = 10;
%viewvec = [1 1 1];
viewvec = [0 0 1];
[X,T] = meshgrid(x-min(x),t);

subplot(1,2,1)
Z = scale*FIcorr(tidx,xidx)/p(1);
surf(X,T,Z)
axis equal;
set(gca, 'DataAspectRatio', [1 16 1]);
xlabel('position (micron)');
ylabel('time (seconds)');
view(viewvec)
title('data');
% to get the same colorscale
Zmin = min(Z(:));
Zmax = max(Z(:));
colormap(jet(256));

ax2 = subplot(1,2,2);
Z = scale*f(p,x,t)/p(1);
surf(X,T,Z)
axis equal;
set(gca, 'DataAspectRatio', [1 16 1]);
xlabel('position (micron)');
ylabel('time (seconds)');
view(viewvec)
title('plot');
set(ax2,'CLim',[Zmin Zmax])

%%
%----------------------------------------
% include initial steepness 
% EXTRACTED D depends strongly on frapTmax
% error bars on D?
%----------------------------------------

lowerIdx = min([FRAPrange notFRAPrange]);
upperIdx = max([FRAPrange notFRAPrange]);

xidx = lowerIdx:upperIdx;
x = xidx*xyresolution;
L = length(x);

frapTmax = 50;%tmax;

tbleach = 1;
tidx = 2:frapTmax;
t = (tidx - tbleach)/framerate;

% function f to be fitted
%
% f = U0*(1 - b*exp(-gamma t)*[1 + erf(m(x-x0)/sqrt(1+4 D m^2 t))]/2 )
%
% where L = sqrt(4 D t)
%
% p = [U0 b gamma x0 m D]

% expfactor = b*exp(-gamma t)/2
expfactor = @(p,x,t) p(2)*repmat(exp(-p(3)*t)', [1 L])/2;

% erfpart = erf(m(x-x0)/sqrt(1+4 D m^2 t))
erfpart = @(p,x,t) erf(1./sqrt(1 + 4*p(5)^2*p(6)*t)'*(x - p(4))*p(5));

% f = U0*(1 - b*exp(-gamma t)*[1 + erf(m(x-x0)/sqrt(1+4 D m^2 t))]/2 )
f = @(p,x,t) p(1)*(1 - expfactor(p,x,t).*(1 + erfpart(p,x,t)));

FIprofile = FIcorr;

% E: energy functional
E = @(p) sum(sum((FIprofile(tidx,xidx) - f(p,x,t)).^2));

U00 = max(FIcorr(:));
b0 = 1;
gamma0 = 1/300;
x0 = mean(x);
m0 = 1;
D0 = 0;

pinit = [U00 b0 gamma0 x0 m0 D0];

% fit 

options = struct();
options.MaxFunEvals = 3000;
[p,Emin,exitflag,output] = fminsearch(E, pinit, options);

% parameters with proper names
param = struct();
param.U00 = p(1);
param.b0 = p(2);
param.gamma0 = p(3);
param.x0 = p(4);
param.m0 = p(5);
param.L0 = 1/param.m0;
param.D = p(6); % compare 0.2 for free Ds, this value is still large
param

allParam.D = param.D;
allParam.L0 = param.L0;
allParam.gamma0 = param.gamma0;

% visualize result

scale = 10;
%viewvec = [1 1 1];
viewvec = [0 0 1];
[X,T] = meshgrid(x-min(x),t);

subplot(1,2,1)
Z = scale*FIprofile(tidx,xidx)/p(1);
surf(X,T,Z)
axis equal;
set(gca, 'DataAspectRatio', [1 16 1]);
xlabel('position (micron)');
ylabel('time (seconds)');
view(viewvec)
title('data');
% to get the same colorscale
Zmin = min(Z(:));
Zmax = max(Z(:));
colormap(jet(256));

ax2 = subplot(1,2,2);
Z = scale*f(p,x,t)/p(1);
surf(X,T,Z)
axis equal;
set(gca, 'DataAspectRatio', [1 16 1]);
xlabel('position (micron)');
ylabel('time (seconds)');
view(viewvec)
title('plot');
set(ax2,'CLim',[Zmin Zmax])


% SAVE RESULTS
save([barefname '_fittedParam'],'allParam');

%% make diffusion constant scatter

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
set(gca,'FontSize', fsize)
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

%%
%----------------------------------------
% initial steepness time for time
%----------------------------------------

% function f to be fitted
%
% f = U0*(1 - A*[1 + erf(m(x-x0)/sqrt(1+4 D m^2 t))]/2 )
%
% p = [U0 A x0 m D]
erfpart = @(p,x,t) erf(1./sqrt(1 + 4*p(4)^2*p(5)*t)'*(x - p(3))*p(4));
f = @(p,x,t) p(1)*(1 - p(2)*(1 + erfpart(p,x,t))/2);

param = struct('U0', [], 'A', [], 'x', [], 'm', [], 'D', []);

for t = 2:tmax

    % E: energy functional
    E = @(p) sum((FI(t,x) - f(p,x,t)).^2);

    U00 = max(FI(t,:));
    A0 = 1;
    x0 = mean(x);
    m0 = 1/5;
    D0 = 0;

     pinit = [U00 A0 x0 m0 D0];

    [p,fminres] = fminsearch(E, pinit);

    param(t).U0  = p(1);
    param(t).A  = p(2);
    param(t).x  = p(3);
    param(t).m = p(4);
    param(t).D = p(5);
end
    
plot(x, FI(t,x))
hold on
N=5;
plot(x, f(p,x), 'r')
hold off

%%
% D is small, 10^(-3) pixel^2/second
% resolution 0.266 mum/pixel 
% so D ~< 10^(-4) mum^2/sec
% compared to 10^-1 - 10^-2 for free Fat

plot([param.D])
mean([param.D])
std([param.D])

%% WHERE AM I

% in the spacetime fit, the fact that the bright side becomes dimmer, even
% after bleach correction, makes the fit want high diffusion, because it is
% the only way to lower the high side, even though the actual edge
% steepness doesn't change
% i should either rescale the steep side more, or find a better data set or
% make something of a time for time fit, but this is similar to rescaling
% the steep side
%
% we could rescale Fat by Ds
%
% also, it is time to start thinking about error bars
%
% on the other hand David is pretty happy with just fitting the width time
% point for time time point



