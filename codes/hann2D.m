function lpFilter = hann2D(fDim, cutoffScale)

if nargin < 2
    cScalex = 0.5;  % cutoff freq in terms of % of kspace FOV
    cScaley = 0.5;
end
if (isscalar(cutoffScale))
    cScalex = cutoffScale;
    cScaley = cutoffScale;
else
    cScalex = cutoffScale(1);
    cScaley = cutoffScale(2);
end
Nx = fDim(1);
Ny = fDim(2);
hannNx = hann(floor(cScalex*Nx));
hannNy = hann(floor(cScaley*Ny));
[tmpX, tmpY] = meshgrid(hannNy, hannNx);
trimHann = tmpX.*tmpY;
[lwx, lwy] = size(trimHann);
lpFilter = zeros(Nx, Ny);
fStartX = floor(Nx/2) + 1 - (floor(lwx/2) + 1);
fStartY = floor(Ny/2) + 1 - (floor(lwy/2) + 1);
lpFilter(fStartX + 1:(fStartX + lwx), fStartY + 1:(fStartY + lwy)) = trimHann;