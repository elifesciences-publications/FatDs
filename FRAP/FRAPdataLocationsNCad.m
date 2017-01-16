datapath = '/Users/idse/Dropbox/Sprinzak/shared/FRAP movies/N-Cad';

movies = {  ...%fullfile('23.6.15','3','3.tif'),...
            fullfile('23.6.15','4','4.tif'),...%fullfile('24.6.15','1','1.tif'),...
            fullfile('24.6.15','2','2.tif'),... 
            fullfile('24.6.15','3','3.tif'),... 
            fullfile('25.6.15','1','1.tif'),...
            fullfile('25.6.15','2','2.tif'),... %5
            fullfile('25.6.15','3','3.tif'),... %6
            fullfile('25.6.15','4','4.tif'),... %7
            fullfile('new movies','2','2.tif'),... %8
            fullfile('new movies','3','3.tif'),... %9
            fullfile('new movies','5','5.tif'),... % 10
            fullfile('new movies','7','7.tif'),... % DONT DO THE BLOCK FOR JUMPS FOR THIS ONE
            fullfile('new movies','13','13.tif'),... % 12 
            fullfile('19.1.16','4','4.tif'),... % 13 
            fullfile('19.1.16','14','14.tif'),...%14
            fullfile('14.2.16','3','3.tif'),... %15
            fullfile('14.2.16','4','4.tif'),... %%fullfile('14.2.16','4','forFigure','movieForFigure8bit.tif'),...%
            fullfile('14.2.16','6','6.tif'),...
            fullfile('14.2.16','8','8.tif'),...
        };

snakes = {  ...%'snakes3 z3.txt',...
            'snakes4 z2.txt',...'snakes1 z3.txt',...
            'snakes2 z2.txt',...
            'snakes3 z2.txt',...
            'snakes1 z2.txt',...
            'snakes2 z2.txt',...
            'snakes3 z1.txt',...
            'snakes4 z2.txt',...
            'snakes.txt',...
            'snakes.txt',...
            'snakes.txt',...
            'snakes.txt',...
            'snakes.txt',...
            'snakes.txt',...
            'snakes.txt',...
            'snakes.txt',...
            'snakes.txt',...
            'snakes.txt',...
            'snakes.txt',...
        };

zslices = {...%3,
            2,....%
            3,2,2,2,2,1,...%,2, 
            2,2,2,1,2,...
            2,2,...
            3,2,3,1};

xyres = {0.27,0.27,0.27,0.27,0.27,0.27,0.27,...
            0.27,0.27,0.27,0.27,0.27,...
            0.27,0.27,0.27,0.27,0.27,0.27};
                                                 
framerates = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,...
                    0.1,0.1,0.1,0.1,0.1,...
                    0.1,0.1,0.1,0.1,0.1,0.1};

dminMax = {NaN,NaN,NaN,NaN,NaN,NaN,NaN,...
                    NaN,NaN,NaN,NaN,NaN,...
                    NaN,NaN,NaN,3,4,5};

cutofftime = {82,NaN,NaN,15,NaN,NaN,NaN,...
                    NaN,30,40,70,65,...
                    48,21,80,NaN,43,68};

windowwidth = {[],[],30,40,[],[],[],...
                        [],45,35,30,30,...
                        [],45,[],[],[],[]};

sign = {NaN,NaN,NaN,NaN,NaN,NaN,NaN,...
                1,-1,-1,1,-1,...
                1,1,...
                1,-1,-1,-1};
            
windowcenter = {24,[],105,40,[],[],[],...
                  [],53,30,65,36,...
                  [],[],[],33,50,46};