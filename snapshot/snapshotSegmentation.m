%--------------------------------------------------------------------------
%
%   snapshot segmentation with stained nuclei
%
%--------------------------------------------------------------------------

clear all;
close all;

% add subfolders of folder containing this script to path
%---------------------------------------------------------
[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

cd(testScriptPath);
cd('../..');
addpath(genpath(pwd));

% import parameters (data and code locations, ...)
%snapshotSegParametersLocal_exp280715;
%snapshotSegParametersLocal_exp020715;
snapshotSegParametersLocal_exp280715_23052017_3

%filepath = {'/Users/idse/Dropbox/Sprinzak/shared/snapshots 07.05.15/6h dox ilastic/dox 1h_5/'};
%corder = {[1 3 2 4]};

% check corder parameter
if numel(corder) ~= numel(filepath)
    warning('corder doesnt have enough elements, assuming default order of channels: 4 1 2 3');
    corder = {};
    for i = 1:numel(filepath), corder{i} = [4 1 2 3]; end
end

[fnames, lims, flabel] = postIlastik(datapath, dataset); 

saveIntermediates = true;

%%
for fi = 1:numel(filepath) %51

    % location of segmentation files
    nucleisegFile   = fullfile(filepath{fi}, 'ilastik', [flabel{fi} '_nuclei_seg.tif']);
    FDFile     = fullfile(filepath{fi}, 'ilastik', [flabel{fi} '_FatDs_seg.tif']);
    bdryFile     = fullfile(filepath{fi}, 'ilastik', [flabel{fi} '_FatDs_bdryseg.tif']);
    
%     if ~exist(nucleisegFile,'file') || ~exist(FDFile,'file')%...
%                                   %  || ~exist(bdryFile,'file')
%         warning(['some segmentation file is missing for ' filepath{fi}]);
%     else
        
    %-----------------------------
    disp('read the data')
    disp(filepath{fi});
    %-----------------------------

    % raw data
    fname = fnames{fi};
    
    % absolute values really don't matter for DIC
    DIC     = mat2gray(imread(fname, corder{fi}(1))); 
    Ds      = imread(fname, corder{fi}(2));
    Fat     = imread(fname, corder{fi}(3));
    Nuc     = imread(fname, corder{fi}(4));
    
    % nuclei segmentation
    nucleiSeg   = imread(nucleisegFile) == 2; % 2 for 28.7.15, 1 for 2.7.15

    % Fat-Ds segmentation
    FDseg     = imread(FDFile);
    DsSeg = FDseg == 1;
    FatSeg = FDseg == 2;
    bg = FDseg == 3;

    % bdry segmentation
    bdryseg     = imread(bdryFile);
    bdries = bdryseg == 2;  % THIS LINE MAY HAVE TO BE EDITED TO 1 

%     % TEMPORARY HACK
%     bdries = FDseg == 3;
%     DsSeg = DsSeg & ~bdries;
%     FatSeg = FatSeg & ~bdries;
%     bg = FDseg == 4 & ~bdries;
%     %-----END HACK------
    
    NSteps = 1;

    % check Ds and Fat are red and green, respectively
    %--------------------------------------------------

    im = cat(3, imadjust(mat2gray(Ds)),imadjust(mat2gray(Fat)),...
                        imadjust(mat2gray(Nuc)));
    %imshow(rawRGB);

    if saveIntermediates
        fname = [num2str(NSteps, '%.2u') '_rawRGB.tif'];
        fullfname = fullfile(filepath{fi},segResultsDir, fname);
        imwrite(im, fullfname);
        NSteps = NSteps + 1;
    end

    % check that segmentation was read correctly
    %--------------------------------------------------

    im = cat(3, DsSeg + bdries, FatSeg + bdries, bg);
    %imshow(bdrySegRGB);

    if saveIntermediates
        fname = [num2str(NSteps, '%.2u') '_bdrysegRGB.tif'];
        fullfname = fullfile(filepath{fi},segResultsDir, fname);
        imwrite(im, fullfname);
        NSteps = NSteps + 1;
    end

    % check segmentation nuclei vs bg
    %--------------------------------------------------

    % nuclei should be yellow on black as much as possible

    im = cat(3, mat2gray(nucleiSeg), nucleiSeg, bg);
    %imshow(nucVsBg);

    if saveIntermediates
        fname = [num2str(NSteps, '%.2u') '_nuclei_bg.tif'];
        fullfname = fullfile(filepath{fi},segResultsDir, fname);
        imwrite(im, fullfname);
        NSteps = NSteps + 1;
    end

    % check segmentation nuclei 
    %--------------------------------------------------

    im = cat(3, mat2gray(nucleiSeg), imadjust(mat2gray(Nuc)), 0*bg);

    if saveIntermediates
        fname = [num2str(NSteps, '%.2u') '_check_nuclei_seg.tif'];
        fullfname = fullfile(filepath{fi},segResultsDir, fname);
        imwrite(im, fullfname);
        NSteps = NSteps + 1;
    end

    %------------------------------------------------------------
    disp('clean up ilastik segmentations (potentially crop)')
    %------------------------------------------------------------

    % crop
    ymax = size(Nuc,1);
    xmax = size(Nuc,2);

    bdriesroi = bdries; %(1:ymax,1:xmax);
    nucroi = nucleiSeg; %(1:ymax,1:xmax);
    bgroi = bg;         %(1:ymax,1:xmax);

    % clean up nuclei segmentation
    nucroi = imerode(nucroi,strel('disk',2));
    nucroi = imopen(nucroi,strel('disk',5));
    nucroi = imfill(nucroi,'holes');

    % old alternative    
    %     nucroi = imclose(nucroi,strel('disk',3));
    %     nucroi = imfill(nucroi, 'holes');
    %     nucroi = bwmorph(nucroi, 'shrink', 3);
    %     nucroi = bwareaopen(nucroi, 50);

    im = cat(3, mat2gray(nucroi), imadjust(mat2gray(Nuc)), 0*nucroi);

    if saveIntermediates
        fname = [num2str(NSteps, '%.2u') '_check_nuclei_seg_clean.tif'];
        fullfname = fullfile(filepath{fi},segResultsDir, fname);
        imwrite(im, fullfname);
        NSteps = NSteps + 1;
    end

    % clean up background
    bgroi = bgroi & ~imdilate(nucleiSeg(1:ymax,1:xmax), strel('disk', 20));
    bgroi = imerode(bgroi, strel('disk', 3));
    bgroi = bwareaopen(bgroi, 100);

    % cleanup boundaries
    cleanBdry = bwareaopen(bdriesroi & ~nucroi, 10);
    cleanBdry = bwmorph(cleanBdry,'thin', 'inf');
    cleanBdry = bwmorph(cleanBdry,'dilate', 1);
    
    %------------------------------------------------------------
    disp('watershedding with nuclei as seeds')
    %------------------------------------------------------------
    
    % impose minima either in background or nuclei
    %Vwater = uint8(bdriesroi);
    %im = imimposemin(Vwater, bgroi | nucroi); %cell
    %im = bwareaopen(smoothDIC < 0.5,100);
    
    % the distance transform provides a smooth potential, where cells get
    % more regular shapes, in particular no 4-connected strings that mess
    % up everything afterwards
    
    im = bwdist(bgroi | nucroi);
    L = watershed(im);

    % at this point there are more regions in L numel(unique(L(:)))-1 
    % than nuclei because of the background regions 

    % give disconnected background regions the same label
    nucLabels = unique(L(nucroi>0));
    bglabels = setdiff(unique(L(:)), [0 nucLabels']);

    for i = 1:length(bglabels)
        L(L==bglabels(i)) = bglabels(1);
    end

    % remove edges between different background regions
    bgedges = imopen(L~=bglabels(1), strel('disk',1));
    idx  = (L~=bglabels(1)) - bgedges > 0;
    newL = L;
    newL(idx) = bglabels(1);

    % at this point #region in newL numel(unique(newL(:)))-1 
    % should be number of nuclei + 1, bc membrane = 0, bg = 1, ..

    finalseg = double(newL==0) + 0.3*double(newL==1);

    % superimpose segmentation of cells and nuclei with Fat/Ds data
    colorsegFatDs = cat(3,  imadjust(mat2gray(Ds)) + finalseg,...
                            imadjust(mat2gray(Fat)) + finalseg,...
                            nucroi + finalseg);

    %figure, imshow(colorsegFatDs);

    im = colorsegFatDs;
    if saveIntermediates
        fname = [num2str(NSteps, '%.2u') '_cells_seg.tif'];
        fullfname = fullfile(filepath{fi},segResultsDir, fname);
        imwrite(im, fullfname);
        NSteps = NSteps + 1;
    end

    % superimpose boundary and cell segementation
    % imshow(cat(3, finalseg + cleanBdry, cleanBdry, nucroi))

    %---------------------------------------------------------------
    disp('remove giant cell cluster messes')
    %---------------------------------------------------------------

    nucCC = bwconncomp(nucroi,8);
    stats = regionprops(nucCC, 'area', 'solidity');
    badNuc = find([stats.Solidity] < 0.9 & [stats.Area] > 2000);
    badCell = uint16(0*badNuc);

    cellsL = newL;
    cellsLclean = cellsL;
    nucroiclean = nucroi;
    
    for i = 1:numel(badNuc)
        % label of bad cells
        badCell(i) = cellsL(nucCC.PixelIdxList{badNuc(i)}(1));
        cellsLclean(cellsL == badCell(i)) = bglabels(1);
        nucroiclean(nucCC.PixelIdxList{badNuc(i)}) = false;
    end
    
    % remove interfaces between bg regions
    newbg = imclose(cellsLclean==1,strel('disk',1));
    cellsLclean(newbg) = 1;
    
    cellsL = cellsLclean;
    nucroi = nucroiclean;

    %---------------------------------------------------------------
    disp('create CellLayer object organizing cell layer properties')
    %---------------------------------------------------------------

    clear cellLayer;
    % as input for a layer that has background, cell and membrane,
    % take image that has membrane 1, cells 0 and background -1
    memseg = int8(cellsL);
    memseg(cellsL == bglabels(1)) = -1;
    memseg(cellsL>1) = 0;
    memseg(cellsL==0) = 1;

    % create a CellLayer object called cellLayer
    cellStateVars = {'FatRaw', 'DsRaw', 'FatSegNuc','DsSegNuc', 'Fat', 'Ds', 'FatExclB','DsExclB'};
    bondStateVars = {'FD','FF','DD','bdry', 'BFI', 'BDI'};
    nTimePts = 1;
    cellLayer = CellLayer(nTimePts, cellStateVars, bondStateVars);

    % initialize timepoint from Lattmin
    % FOR SOME REASON minVertexDist == 0 WILL MESS UP BONDLABELS
    % I SHOULD FIGURE OUT WHY
    t = 1;
    options = struct('closeSize', 0, 'minVertexDist', 1, 'bgCells', 0, 'trim', 0);
    cellLayer.initTime(t, 'image', memseg, options);
    
    %-----------------------------------------------
    disp('determine Fat/Ds levels in each cell')
    %-----------------------------------------------

    cellL = cellLayer.cellL(t);

    %bdrypixel = find(bdriesroi);
    CC = bwconncomp(cellL>0,4);
	nucCC = bwconncomp(nucroi,8);

    FatNoB = double(Fat);
    DsNoB = double(Ds);
    FatNoB(bdriesroi) = NaN;
    DsNoB(bdriesroi) = NaN;
    
    if nucCC.NumObjects~=CC.NumObjects
        warning('Ncells ~= Nnuclei!');
    else

    for i = 1:CC.NumObjects

        % raw intensity in cell mask
        FatRaw = mean(Fat(CC.PixelIdxList{i}));
        DsRaw = mean(Ds(CC.PixelIdxList{i}));
        
        FatExclB = nanmean(FatNoB(CC.PixelIdxList{i}));
        DsExclB = nanmean(DsNoB(CC.PixelIdxList{i}));
        
%         % slow but freedom to change the cellmask
%         mask = false(size(FatSeg));
%         mask(CC.PixelIdxList{i}) = true;
%         mask = imerode(mask,strel('disk',4));
%         FatRaw = mean(Fat(mask));
%         DsRaw = mean(Ds(mask));

        % update cellLayer
        s = cellLayer.cells{t}(i).state;
        s(1:2) = [FatRaw DsRaw];
        s(7:8) = [FatExclB DsExclB];
        cellLayer.cells{t}(i).setState(s);
        
        % segmented amount of red/green in nuclear mask
        FatSegNuc =  mean(FatSeg(nucCC.PixelIdxList{i}));
        DsSegNuc = mean(DsSeg(nucCC.PixelIdxList{i}));

        % there is no guarantee that nuclear and cell labels match even
        % thought their total number is the same
        nucCell = cellLayer.L(nucCC.PixelIdxList{i}(1));

        % update cellLayer
        s = cellLayer.cells{t}(nucCell).state;
        s(3:4) = [FatSegNuc DsSegNuc];
        cellLayer.cells{t}(nucCell).setState(s);
    end

    %-----------------------------------------------
    disp('threshold segmented nuclear levels')
    %-----------------------------------------------

    cellsmask = imerode(newL>1,strel('disk',1));
    FatCellsMask = false(size(cellsmask));
    DsCellsMask = false(size(cellsmask));

    cutoff = 0;
    tic
    for i = 1:CC.NumObjects

        cstate = cellLayer.cells{t}(i).state;
        FatLevel = cstate(3);
        DsLevel = cstate(4);
        relLevel= (FatLevel-DsLevel)/(FatLevel + DsLevel);

        if relLevel > cutoff

            cstate(5:6) = [1 0];
            FatCellsMask(CC.PixelIdxList{i}) = true;

        elseif relLevel < -cutoff 

            cstate(5:6) = [0 1];
            DsCellsMask(CC.PixelIdxList{i}) = true;

        else
            cstate(5:6) = [0 0];
        end

        cellLayer.getCell(t,i).setState(cstate);
    end
    FatCellsMask = FatCellsMask & ~nucroi;
    DsCellsMask = DsCellsMask & ~nucroi;
    otherCellsMask = cellsmask & ~FatCellsMask &~ DsCellsMask;
    toc

    segmented = cat(3, mat2gray(DsCellsMask) + cleanBdry, FatCellsMask + cleanBdry, otherCellsMask & ~cleanBdry);
%     figure,
%     imshow(segmented);

    im = segmented;
    if saveIntermediates
        fname = [num2str(NSteps, '%.2u') '_cells_colored_nucMask.tif'];
        fullfname = fullfile(filepath{fi},segResultsDir, fname);
        imwrite(im, fullfname);
        NSteps = NSteps + 1;
    end
    % %%
    % imshow(segmented);
    % options = struct('cellIndex', false, 'transparent', true, 'lineWidth', 2,...
    %                     'edgeColor','m');
    % cellLayer.visualize(t, options)

    %--------------------------------------------------------------
    disp('interface counting')
    %--------------------------------------------------------------

    FatDsIfMask = false(size(cellsmask));

    % for each bond
    for i = 1:numel(cellLayer.bonds{t})

        bondCellIdx = cellLayer.bonds{t}(i).cellInd;

        % if it not an outside bond (if it is a cell-cell interface)
        if ~any(bondCellIdx == 0)

            states = cat(1,cellLayer.cells{t}(bondCellIdx).state);
            states = states(:,5:6);

            % Fat-Ds interface
            if all(sum(states) == [1 1])

                % set CanHaz to true, Haz to NaN until determined
                cellLayer.bonds{t}(i).setState([true false false NaN NaN NaN]);
                blabel = cellLayer.bonds{t}(i).label;
                FatDsIfMask(cellLayer.L == blabel) = true;

            % Ds-Ds
            elseif all(sum(states) == [0 2])
                cellLayer.bonds{t}(i).setState([false false true NaN NaN NaN]);

            % Fat-Fat
            elseif all(sum(states) == [2 0])
                cellLayer.bonds{t}(i).setState([false true false NaN NaN NaN]);
            end
        end
    end

    dilCanHaz = imdilate(FatDsIfMask, strel('disk',3));
    %FatDsMask = imdilate(FatDsMask, strel('disk',1));
    potBdries = cat(3,  mat2gray(DsCellsMask & ~dilCanHaz) + dilCanHaz,...
                        (FatCellsMask & ~dilCanHaz) + dilCanHaz,...
                        (otherCellsMask & ~dilCanHaz)) ;

    %imshow(potBdries);
    im = potBdries;
    if saveIntermediates
        fname = [num2str(NSteps, '%.2u') '_potential_bdries.tif'];
        fullfname = fullfile(filepath{fi},segResultsDir, fname);
        imwrite(im, fullfname);
        NSteps = NSteps + 1;
    end
   
    %--------------------------------------------------------------
    disp('identify Fat-Ds interfaces that actually form boundaries')
    %--------------------------------------------------------------

    % watershed the boundaries
    bwif = cellLayer.L > numel(cellLayer.cells{1});
    im = bwdist(bwif);
    L = watershed(im);
    
    % visualize
%     Lrgb = label2rgb(L);
%     for c=1:3
%         Lrgb(:,:,c) = double(~bwif).*double(Lrgb(:,:,c));
%     end
%     imshow(Lrgb,[])
    
    % construct a lookup table between watershed L and interface L
    ifCC = bwconncomp(bwif);
    LUT = zeros([ifCC.NumObjects 2],'uint16');

    for i = 1:ifCC.NumObjects
        LUT(i,:) = [L(ifCC.PixelIdxList{i}(1)) cellLayer.L(ifCC.PixelIdxList{i}(1))];
    end
    LUT = sort(LUT);

    % to haz or not to haz
    Haz = false(size(bdriesroi));
    
    bdryCC = bwconncomp(cleanBdry);
    
    for bi = 1:bdryCC.NumObjects

        bdryPL = bdryCC.PixelIdxList{bi};
        
        basinvals = L(bdryPL);

        basinvals = basinvals(basinvals>0);
        basinidx = unique(basinvals);
        n = histc(basinvals,basinidx);
        
        % THIS IS A TUNABLE PARAMETER
        % minimum number of accumulation pixels to count it toward some
        % interface
        basinidx = basinidx(n > 20);
        
%         R = uint8(cellLayer.L);% & cellLayer.L < numel(cellLayer.cells{1}));
%         G = cellLayer.L == bondLabel;
%         B = G;
%         imshow(255*cat(3,R,G,B))
%         

        if ~isempty(basinidx)
            
            for i = 1:numel(basinidx)
                
                bondLabel = LUT(basinidx(i),2);
                opposingBonds = cellLayer.getBond(t,bondLabel);
            
                % ~isempty: temporary workaround until I deal with isolated pairs
                % state(1) == 1: only assign accumulation to interfaces between Fat
                % and Ds, not other interfaces
                if ~isempty(opposingBonds) && opposingBonds(1).state(1)  == 1  

                    % intensity of boundary
                    BFI = mean(Fat(bdryPL));
                    BDI = mean(Ds(bdryPL));

                    % a two sided cell sharing one interface with another
                    % cell and one with the outside has two opposing bonds 
                    % but only one can be accumulating boundary, the if
                    % statement excludes the one facing the outside (
                    if all(opposingBonds(1).cellInd > 0)
                        opposingBonds(1).setState([true false false true BFI BDI]);
                    end
                    if all(opposingBonds(2).cellInd > 0)
                        opposingBonds(2).setState([true false false true BFI BDI]);
                    end

                    Haz(ifCC.PixelIdxList{basinidx(i)}) = true;
                end
            end
        end
    end

    % image 
    dilHaz = imdilate(Haz, strel('disk',3));
    %FatDsMask = imdilate(FatDsMask, strel('disk',1));
    segmented = cat(3,  mat2gray(DsCellsMask & ~dilHaz) + dilHaz,...
                        (FatCellsMask & ~dilHaz) + dilHaz,...
                        (otherCellsMask & ~dilHaz)) ;
%     imshow(segmented)
% 
%     segmented = cat(3,  mat2gray(DsCellsMask & ~cleanBdry) + cleanBdry,...
%                         (FatCellsMask & ~cleanBdry) + cleanBdry,...
%                         (otherCellsMask & ~cleanBdry)) ;
%     figure, imshow(segmented)
                    
    im = segmented;
    %if saveIntermediates
        fname = '01_seg_full_new.tif';
        fullfname = fullfile(filepath{fi},segResultsDir, fname);
        imwrite(im, fullfname);
        NSteps = NSteps + 1;
    %end

% %%
%     %--------------------------------------------------------------
%     disp('identify Fat-Ds interfaces that actually form boundaries')
%     %--------------------------------------------------------------
% 
%     % to haz or not to haz
%     Haz = false(size(bdriesroi));
% 
%     ifCC = bwconncomp(FatDsIfMask);
%     bdryCC = bwconncomp(bdriesroi);
%     bdryL = bwlabel(bdriesroi);
%     
%     % make a list of assigned boundaries
%     % then we can count those without overcounting double assignement
%     assignedBdry = [];
%     
%     for bi = 1:ifCC.NumObjects
% 
%         intersectSize = sum(bdryL(ifCC.PixelIdxList{bi})>0);
%         bondLabel = cellLayer.L(ifCC.PixelIdxList{bi}(1));
%         
%         % more than one boundary can be intersected
%         % for now i just combine them
%         % there can be overcounting due to this
%         bdryLabel = setdiff(unique(bdryL(ifCC.PixelIdxList{bi})),0);
%         
%         %if intersectSize./bdrySize > 0.2 % THRESHOLD TO PLAY WITH
%         if intersectSize > 0
% 
%             assignedBdry = [assignedBdry(:)' bdryLabel(:)'];
%             opposingBonds = cellLayer.getBond(t,bondLabel);
% 
%             % temporary workaround until I deal with isolated pairs
%             if ~isempty(opposingBonds)
% 
%                 BFI = mean(Fat(cat(1,bdryCC.PixelIdxList{bdryLabel})));
%                 BDI = mean(Ds(cat(1,bdryCC.PixelIdxList{bdryLabel})));
%                 
%                 opposingBonds(1).setState([true false false true BFI BDI]);
%                 opposingBonds(2).setState([true false false true BFI BDI]);
%             end
% 
%             Haz(ifCC.PixelIdxList{bi}) = true;
%         end
%     end
% 
%     % boundary intensities                  
%     assignedBdry = unique(assignedBdry);
%     bdryFI = zeros(size(assignedBdry));
%     bdryDI = zeros(size(assignedBdry));
%     for i = 1:numel(bdryFI)
%         bdryFI(i)=mean(Fat(bdryCC.PixelIdxList{assignedBdry(i)}));
%         bdryDI(i)=mean(Ds(bdryCC.PixelIdxList{assignedBdry(i)}));
%     end
% 
%     % image 
%     dilHaz = imdilate(Haz, strel('disk',3));
%     %FatDsMask = imdilate(FatDsMask, strel('disk',1));
%     segmented = cat(3,  mat2gray(DsCellsMask & ~dilHaz) + dilHaz,...
%                         (FatCellsMask & ~dilHaz) + dilHaz,...
%                         (otherCellsMask & ~dilHaz)) ;
%     imshow(segmented)
% 
%     segmented = cat(3,  mat2gray(DsCellsMask & ~cleanBdry) + cleanBdry,...
%                         (FatCellsMask & ~cleanBdry) + cleanBdry,...
%                         (otherCellsMask & ~cleanBdry)) ;
%     figure, imshow(segmented)
%                     
%     im = segmented;
%     if saveIntermediates
%         fname = [num2str(NSteps, '%.2u') '_seg_full.tif'];
%         fullfname = fullfile(filepath{fi},segResultsDir, fname);
%         imwrite(im, fullfname);
%         NSteps = NSteps + 1;
%     end
    
    %--------------------------------------------------------------
    disp('save results for later analysis')
    %--------------------------------------------------------------

%     save(fullfile(filepath{fi},segResultsDir,[flabel{fi} '_seg617']), 'cellLayer',...
%                                     'bdryFI', 'bdryDI', 'assignedBdry')

    save(fullfile(filepath{fi},segResultsDir,[flabel{fi} '_seg623']), 'cellLayer');

    end
    %end
end
