clear all;
close all;

% add subfolders of folder containing this script to path
[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

cd(testScriptPath);
cd('..');
addpath(genpath(pwd));

dataDir = '/Users/idse/Dropbox/Sprinzak/shared/polarity/160410_beads';
cd(dataDir);

%%

fname = 'f d beads 11.4.16.lif - Series010.tif';
imz = readStack('.',fname);

%imz  = imz(400:900,100:600,:,:);

%% make beads mask

zidx = 18;
lim = stretchlim(mat2gray(imz(:,:,zidx,5)));
beadsmask = mat2gray(imz(:,:,zidx,5)) > lim(2);
beadsmask = imerode(beadsmask, strel('disk',1));
beadsmask = bwareaopen(beadsmask,3);
beadsmask = imclearborder(beadsmask);
beadsmask = imdilate(beadsmask, strel('disk',5));
imshow(beadsmask,[]);

%% overlay on red and green

beadsR = beadsmask.*imadjust(mat2gray(imz(:,:,zidx,1)));
beadsG = beadsmask.*imadjust(mat2gray(imz(:,:,zidx,2)));
beadsRG = cat(3,beadsR,beadsG,0*beadsG);
imshow(beadsRG)

%%

stats = regionprops(beadsmask, 'centroid', 'boundingbox');
CM = cat(1,stats.Centroid);

margin = 4;
scale = 10;
shift = zeros([numel(stats) 2]);
imAfter = {};
imBefore = {};

for i = 1:numel(stats)
    
    b = stats(i).BoundingBox;
    b(1:2) = b(1:2) - margin;
    b(3:4) = b(3:4) + 2*margin;
    ymax = uint32(b(2)+b(4));
    xmax = uint32(b(1)+b(3));
    
    if ymax < size(imz,1) && xmax < size(imz,2)
        
        imR = mat2gray(imz(b(2):ymax,b(1):xmax,zidx,1));
        imG = mat2gray(imz(b(2):ymax,b(1):xmax,zidx,2));
        %imshow(cat(3,imR,imG,0*imG),[])

        [shiftx,shifty] = xcorr2fft(imresize(imR,scale),imresize(imG,scale));
        shiftx = shiftx/scale;
        shifty = shifty/scale;
        shift(i,:) = [shiftx shifty];
        
    else
        shift(i,:) = [NaN NaN];
    end
    
    xform = [1  0  0; 0  1  0; shiftx shifty  1];
    tform = maketform('affine',xform);
    imRp = imtransform(imR,tform,'XData',[1 size(imR,2)],'YData',[1 size(imR,1)]);
    
    imAfter{i} = cat(3,imRp,imG,0*imR);
    imBefore{i} = cat(3,imR,imG,0*imR);
end

CM = CM(~isnan(shift(:,1)),:);
shift = shift(~isnan(shift(:,1)),:);
meanShift = mean(shift);

%% visualize individual bead correction

i = 6
figure(1), imshow(imBefore{i});
figure(2), imshow(imAfter{i});

%% visualize displacement between channels

s = 100;
imshow(cat(3,imadjust(mat2gray(imz(:,:,zidx,1))),imadjust(mat2gray(imz(:,:,zidx,2))),0*imadjust(mat2gray(imz(:,:,zidx,1)))))
hold on
quiver(CM(:,1),CM(:,2),s*shift(:,1),s*shift(:,2),0)
quiver(512,512,s*meanShift(:,1),s*meanShift(:,2),0,'m')
hold off
saveas(gcf,'displacements.png');

%% visualize displacement between channels

xform = [1  0  0; 0  1  0; meanShift(1) meanShift(2)  1];
tform = maketform('affine',xform);
imR = imz(:,:,zidx,1);
imG = imz(:,:,zidx,2);
imRp = imtransform(imR,tform,'XData',[1 size(imR,2)],'YData',[1 size(imR,1)]);

%%
figure, imshow(cat(3,mat2gray(imRp),mat2gray(imG),0*mat2gray(imR)))

imwrite(imG, 'corrected.tif')
imwrite(imRp, 'corrected.tif', 'WriteMode','Append')

imwrite(imG, 'uncorrected.tif')
imwrite(imR, 'uncorrected.tif', 'WriteMode','Append')
