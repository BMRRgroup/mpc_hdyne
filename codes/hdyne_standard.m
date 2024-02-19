function [img] = hdyne_standard(pk, nover)
  %HOMODYNE Reconstruct partial k-space data using the homodyne algorithm
  %
  %  IMG = hdyne_standard(PK,NOVER) reconstructs the real image IMG from the
  %  2D partial k-space data PK where NOVER samples were acquired prior to the
  %  echo.  The parameter NTRAN determines the sharpness of the transitions in
  %  the weighting functions that extract the symmetric and asymmetric parts of
  %  PK.  The echo is assumed to be in the lower half of PK; ie, echoes form at
  %  row N where N=NY-NOVER and NY is the number of rows in PK.
  %
  %  The fractional echo data is assume in the column dimension
  %

  ntran = 2;
  datatype = 'k';
  
  
  ksize = size(pk);
  ksize(end+1:6) = 1;

  if ismatrix(pk)
    ksize(3) = 1;
  end
    
  np = ksize(1);
  nx = ksize(2);
  nf = 2 * (np - nover);
  
  % form the weighting matrices--only 1D, nf in length
  y = 1:nf;
  q = np - 2*ntran;
  p = nf - q;
  u = 2 - 1./(1+exp((p-y)/ntran)) - 1./(1+exp((q-y)/ntran));
  v = 1./(1+exp((y-q)/ntran)) - 1./(1+exp((y-p)/ntran));
  
  pk(nf,:) = 0;  % zeropad the partial k-space data. Works for ndims > 2 too !

  u = repmat(u', 1, nx);  % asymmetric weights
  v = repmat(v', 1, nx);  % symmetric weights
  au = zeros(size(pk));
  av = zeros(size(pk));
  for k = 1:prod(ksize(3:end))
      au(:,:,k) = pk(:,:,k) .* u;  % asymmetric k-space data
      av(:,:,k) = pk(:,:,k) .* v;  % symmetric k-space data
  end
  clear  u v
 
  % get image domain data
  a = fft2d(au,-1);  % asymmetric image data
  s = fft2d(av,-1);  % symmetric image data
  clear au av
  img = real(a .* exp(-1i*angle(s)));  % image is real part of phase-corrected data
end







