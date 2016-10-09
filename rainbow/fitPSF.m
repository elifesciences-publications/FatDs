%----------------------------------------
% fit Gaussian to PSF
%----------------------------------------


% define location of segmented data and read
%--------------------------------------------
dataDir  = '/Users/idse/Dropbox/Sprinzak/shared/images for psf 13.7.15/';

fname = {'1.7.15_Series029_prettygood.tif',...
         '1.7.15_Series028_prettygoodbutartefacts.tif',...
         '1.7.15_Series042_underexposed.tif',...
         '1.7.15_Series022_saturated.tif'};
xyres = {0.0212, 0.0212, 0.0436, 0.0655};
zres = 0.3002;

fi = 1;
info = imfinfo(fullfile(dataDir,fname{fi}));
nSlices = numel(info);
im = zeros([info(1).Height, info(1).Width nSlices]);

for i = 1:nSlices
	im(:,:,i) = imread( fullfile(dataDir,fname{fi}), i);
end

% %%
% % fit 2D gaussian to slice
% %--------------------------------------------
% 
% slice = im(:,:,10);
% %imshow(slice,[])
% 
% A = max(slice(:));
% sig = 6;
% x0 = size(slice,2)/2;
% y0 = size(slice,1)/2;
% 
% pinit = [A sig x0 y0];
% gauss2d  = @(p,x,y) p(1)*exp(-((x-p(3)).^2+(y-p(4)).^2)/(2*p(2)^2));
% 
% [X,Y] = meshgrid(1:size(slice,2),1:size(slice,1));
% 
% imshow(gauss2d(pinit, X,Y),[])
% 
% % E: energy functional
% E= @(p) sum(sum((slice - gauss2d(p,X,Y)).^2));
% 
% % fit 
% options = struct();
% options.MaxFunEvals = 3000;
% [p,Emin,exitflag,output] = fminsearch(E, pinit, options);
% 
% p
% imshow(cat(3,slice/p(1),gauss2d(p,X,Y)/p(1),slice*0))

% fit 3D gaussian
%--------------------------------------------

slice = im(:,:,10);
%imshow(slice,[])

A = max(slice(:));
sig = 6*(0.0212/xyres{fi});
x0 = size(slice,2)/2;
y0 = size(slice,1)/2;
sigz = 2;
z0 = 8;

pinit = [A sig x0 y0 sigz z0];
gauss3d  = @(p,x,y,z) gauss2d(p,x,y).*exp(-(z-p(6)).^2/(2*p(5)^2));

[X,Y,Z] = meshgrid(1:size(im,2),1:size(im,1),1:size(im,3));

%imshow(gauss2d(pinit, X,Y),[])

% E: energy functional
E= @(p) sum(sum(sum((im - gauss3d(p,X,Y,Z)).^2)));

% fit 
options = struct();
options.MaxFunEvals = 3000;
[p,Emin,exitflag,output] = fminsearch(E, pinit, options);

p
PSFfitted = gauss3d(p,X,Y,Z);

x = round(p(3));
z = round(p(5));

imshow(squeeze(PSFfitted(:,x,:)),[])
imshow(squeeze(PSFfitted(:,:,z)),[])

sigz = p(5)*zres
sigxy = p(2)*xyres{fi}

