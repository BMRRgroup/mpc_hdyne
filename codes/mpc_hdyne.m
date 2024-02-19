function [img_pF, diff_img] = mpc_hdyne

% Example codes for the homodyne algorithms from the manuscript entitled 
% "Motion-induced phase corrected homodyne reconstruction for partial Fourier
% single-shot diffusion-weighted echo planar imaging of the liver"
%
% Anh Van
% BMRR Group, Technical Unitersity of Munich
% 02/2024

load ../data/img_fF.mat img_fF;

addLphase = 1;
ntyp = 4;     % no phase, with phase, with phase + shift to full, with phase + shift to crop

[Ny, Nx] = size(img_fF);

svol = zeros(Ny, Nx, ntyp);
svol(:,:, 1) = abs(img_fF);     % no phase simulation
svol(:,:, 2) = img_fF;          % with native motion-induced phase (MiPe)
                        
if (addLphase)
    % apply additional linear phase
    [~, ygrid] = meshgrid((1:Nx), (1:Ny));
    full_slopey = 2*pi*14/Ny;    % + shift to full 
    crop_slopey = -2*pi*10/Ny;   % - shift to crop
        
    lphase = exp(1j*full_slopey*ygrid);
    svol(:, :, 3) = svol(:, :, 2).*lphase;  % MiPe + shift-to-full
       
    lphase = exp(1j*crop_slopey*ygrid);
    svol(:, :, 4) = svol(:, :, 2).*lphase;  % MiPe + shift-to-crop
end

kvol = fft2d(svol, 1);
kRange = [-Ny/2; Ny/2 - 1];

pFyall = 5/8;  %[5/8; 6/8; 7/8];        
npF = length(pFyall);
npc = 3;   %  homodyne, adaptive homodyne, mpc-hdyne

% initialize variables
img_pF = zeros(Ny, Nx, ntyp, npc, npF);   % reconstructed pF images
diff_img = zeros(size(img_pF));           % difference images with full Fourier

for ipF = 1:npF
    pFy = pFyall(ipF);    
    pF = [1; pFy];
    
    % simulate pFourier acquisition    
    [pkvol, pFkRange, ~] = simpFourier(permute(kvol, [2 1 3 4]), kRange, pF, 1);
    pkvol = permute(pkvol, [2 1 3 4]);   % zero-filled pF data
    
    %....... Standard homodyne......
    hnover = - pFkRange(1);
    pFStart = Ny/2 + 1 - hnover;
    tmpk = pkvol(pFStart:end, :, :, :);
    tmpk = flip(tmpk, 1);
    img_pF(:, :, :, 1, ipF) = hdyne_standard(tmpk, hnover); 
    
    clear tmpk;
    %........Adaptive homodyne and POCS...... 
     % - Linear phase is corrected by shifting k-space and adjust k-space
     % extent (partial Fourier factor). 
     % - Non-linear phase is correct by phase subtraction
     kcenter = floor(Ny/2) + 1;  % pkvol is zero-filled pF data
     Nky = Ny;
     for ityp = 1:ntyp
         
             tmpk = pkvol(:, :, ityp);   
             
             %.....find k-space center and shift k-space
             [cy , ~] = find(abs(tmpk) == max(abs(tmpk(:))));
             cy = cy - kcenter;  
             if (cy < 0)    % kcenter shifted to crop region, need to extend k-space when shifting back
                 newk = zerofill(tmpk, [Nky - 2*cy, size(tmpk, 2)], [0, 0]);
             else
                 newk = tmpk;
             end
             newk = circshift(newk, -cy, 1);

             if (cy < 0)
                 % FIXME: truncate back to original k-space size. Can we do better?
                 newk = newk((-cy + 1): (Nky - cy), :); 
             end

             %......adjusted partial Fourier parameters
             hnover = - pFkRange(1) + cy;
             pFStart = Ny/2 + 1 - hnover; 

             %......adaptive homodyne.............
             pFkvol = newk(pFStart:end, :);
             pFkvol = flip(pFkvol, 1);
             img_pF(:, :, ityp, 2, ipF) = hdyne_standard(pFkvol, hnover); 


             
             %.....Non-linear phase correction for homodyne
             shiftedImg = fft2d(newk, -1);               
             cimg = remove_lrPhase(shiftedImg, (-pFkRange(1) + cy)/Ny);
             ckvol = fft2d(cimg, 1);
             
             clear tmpk newk shiftedImg;            
             
             %....for mpc-hdyne.........
             pFkvol = ckvol(pFStart:end, :);
             pFkvol = flip(pFkvol, 1);
             img_pF(:, :, ityp, 3, ipF) = hdyne_standard(pFkvol, hnover);              
             
             clear ckvol pFkvol;
         
     end     
end
img_pF = flip(abs(img_pF), 1);
img_pF = circshift(img_pF, 1, 1);

% difference with full Fourier
for ipF = 1:npF
    for ipc = 1:npc
        for ityp = 1:ntyp
            diff_img(:, :, ityp, ipc, ipF) = abs(img_pF(:, :, ityp, ipc, ipF)) - abs(img_fF);
        end
    end
end