% data and code locations for snapshot segmentation
%--------------------------------------------------

%epitheliumRepo = 'D:\My Documents\My Desktop\Idse_2015\epitheliumAnalysis';
epitheliumRepo = '/Users/idse/repos/epitheliumAnalysis';

% datapath = 'D:\Olga\DoxVsNoDoxJun24';
% %filepath{fi} = {fullfile(datapath, 'dox 6h_9')};
%  filepath = {fullfile(datapath, '1_12h_0'),...
%              fullfile(datapath, '1_12h_1'),...
%              fullfile(datapath, '1_12h_2'),...
%              fullfile(datapath, '1_12h_3'),...
%              fullfile(datapath, '1_nodox_1'),...
%              fullfile(datapath, '1_nodox_2'),... 
%              fullfile(datapath, '1_nodox_3')};

datapath = '/Users/idse/Dropbox/Sprinzak/data';
filepath = {fullfile(datapath, 'dox on', 'dox 6h', '2_10'),...
            fullfile(datapath, 'no dox', 'no dox_10')};
        
% colors R=1, G=2, etc, channels: [DIC Ds Fat Nuc]
% so [4 1 2 3] makes Ds=R, Fat=G, Nuc=B
% corder = {[4 1 2 3], [4 1 2 3]};
        
corder = {};%[1 2 3 4],[1 2 3 4],[1 2 3 4],[4 2 1 3]};

% save intermediate results for quality inspection
%--------------------------------------------------
saveIntermediates = true;
segResultsDir = 'matlab_seg';