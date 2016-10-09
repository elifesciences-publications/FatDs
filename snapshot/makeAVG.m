function makeAVG(dataDir,fname,z1,z2)

    info = imfinfo(fullfile(dataDir,fname));
    nSlices = numel(info);
    im = zeros([info(1).Height, info(1).Width 10 4],'uint16');

    for i = 1:nSlices
        z = floor(i/4) + 1;
        c = rem(i,4);
        if c==0, c=4; end
        im(:,:,z,c) = imread( fullfile(dataDir,fname), i);
    end

    avgim = uint16(squeeze(mean(im(:,:,z1:z2,:),3)));
    imwrite(avgim(:,:,1), fullfile(dataDir,['AVG_' fname]));
    for i = 2:4
        imwrite(avgim(:,:,i), fullfile(dataDir,['AVG_' fname]),'WriteMode','append');
    end
end