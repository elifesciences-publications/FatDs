function dX = register1D(f1, f2, upsample)
    % register 1D datasets, with subpixel accuracy
    % 
    % dX = register1D(f1, f2)
    % dX = register1D(f1, f2, upsampling)
    %
    % f1, f2:       m by 1 matrices 
    % upsampling:   upsampling factor in fft for subpixel accuracy
    %               an odd integer is better
    % dX:           shift between f2 and f1
    %
    % Idse Heemskerk, Aug 2014
    
    if nargin < 3
        upsample = 1;
    end
    
    % Fourier tranform of the cross-correlation
    F = fftshift(fft(f1));
    Fc = fftshift(conj(fft(f2)));
    R0 = F.*Fc;
    
    % for subpixel accuracy, zero pad the Fourier transform of the
    % correlation
    L = size(R0,1);
    R = zeros([upsample*L 1]);
    lowerIdx = round((upsample-1)*L/2 + 1);
    
    R(lowerIdx:lowerIdx + L-1) = R0;

    
    % get the shift
    c = ifftn(ifftshift(R));
    [~,idx] = max(c(:));
    
    dX = (idx - 1)/upsample;
    
    if L - dX < dX
        dX = dX - L;
    end
end