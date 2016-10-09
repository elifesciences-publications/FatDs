% writeTiff(I, filename)
% writeTiff(R,G,B filename)

function writeTiff(varargin)

if nargin==2 
    imagestack = varargin{1};
    barefilename = varargin{2};
elseif nargin==4
    imagestack(:,:,:,1) = varargin{1};
    imagestack(:,:,:,2) = varargin{2};
    imagestack(:,:,:,3) = varargin{3};
    barefilename = varargin{4};
end

%imwrite only writes doubles to tiff but don't convert the whole stack to
%tiff at once because that would be 10^30x8 = 8 GB!!

% Compression none is essential for Fiji to be able to read these
imwrite(double(squeeze(imagestack(:,:,1,:))), strcat(barefilename,'.tif'),'Compression','none');
for j=2:size(imagestack,3)
   imwrite(double(squeeze(imagestack(:,:,j,:))), strcat(barefilename,'.tif'),'Compression','none', 'writemode', 'append');
end
