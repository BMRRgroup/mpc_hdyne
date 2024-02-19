function img_out = fft2d(img_in, opt)
% opt = 1: perform 2D FFT on input with dimension >= 2
% opt = -1: perform 2D iFFT on input with dimension >= 2
if (opt == 1)
    img_out = fftshift(fft(ifftshift(img_in, 1),[], 1), 1);
    img_out = fftshift(fft(ifftshift(img_out, 2),[], 2), 2);
else
    if (opt == -1)
        img_out = fftshift(ifft(ifftshift(img_in, 1),[], 1), 1);
        img_out = fftshift(ifft(ifftshift(img_out, 2),[], 2), 2);
    end
end
