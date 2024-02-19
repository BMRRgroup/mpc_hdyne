function out_data = remove_lrPhase(in_data, r_lrFact)
% remove the low resolution phase from input data
% the low resolution region is defined by r_lrFact

[Nx, Ny, nvol] = size(in_data);
tmp = fft2d(in_data, 1);

% get low res phase
lpFilter = hann2D([Nx, Ny], [0.9, 2*r_lrFact]);
lpFilter = repmat(lpFilter, [1, 1, nvol]);
lr_in_data = lpFilter.*tmp;
lr_in_data = fft2d(lr_in_data, -1);
% lr_in_data = lr_in_data./abs(lr_in_data);
% lr_in_data(isinf(lr_in_data)) = 0;
% lr_in_data(isnan(lr_in_data)) = 0;

%out_data = in_data.*conj(lr_in_data);

out_data = in_data.*exp(-1j*angle(lr_in_data));


