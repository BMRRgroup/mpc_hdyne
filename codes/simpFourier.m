function [pFkdat, pFkRange, pF] = simpFourier(kdat, kRange, pF, returnzfill)

% Simulate partial Fourier from full Fourier data in the 2nd dimension

if nargin < 4
    returnzfill = 1; % return zerofilled k-space
end

pFdir = find(pF < 1);
pFscalar = pF(pFdir);
nKfull = kRange(2) - kRange(1) + 1;          % full kspace size
kOver = floor(nKfull*pFscalar - kRange(2));  % number of over scan
pFscalar = (kOver + kRange(2) + 1)/nKfull;   % adjusted partial Fourier factor



kSize = size(kdat);
kdat = reshape(kdat, kSize(1), kSize(2), []);
pFkdat = kdat;
if (pFdir == 2)   % for now only support ky partial fourier simulation
    pF = [1 pFscalar 1];
    pFkRange = [-kOver kRange(2)];
    if (returnzfill)
        pFkdat(:, (1:abs(kRange(2)) - kOver), :) = 0;
    else
        tmpk = kdat(:, (abs(kRange(1)) - kOver + 1):end, :);
        pFkdat = tmpk;
        kSize(2) = size(pFkdat, 2);
    end
end

% with filter for ringing control
if (returnzfill)
    nf = nKfull;
    np = floor(nKfull/2) + kOver;
    ntran = 1;
    y = 1:nf;
    q = np - 2*ntran;
    p = nf - 1;
    v = 1./(1+exp((y-q)/ntran));
    v = flip(v, 2);
    v = repmat(v, size(kdat, 1), 1, size(pFkdat, 3));
    pFkdat = pFkdat.*v;
    clear v;
end

pFkdat = reshape(pFkdat, kSize);