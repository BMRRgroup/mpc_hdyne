function b = zerofill(a, dim, loc)
%ZEROFILL Extend the size of a 2D matrix by adding zeros.
%
%  B = ZEROFILL(A, DIM) returns the 2D matrix A zerofilled so that
%  SIZE(B)=DIM.  Zeros are added symmetrically to the top and bottom
%  of A as needed so that the output matrix has the specified size (DIM).
%  If 0 is given for one or both elements of DIM, the respective dimension
%  is not zerofilled.
%
%  B = ZEROFILL(A, DIM, LOC) controls where the zerofilling occurs:
%    if LOC(1) = -1, DIM(1) zeros are added at the top;
%    if LOC(1) = 0, DIM(1)/2 zeros are added at the top and bottom;
%    if LOC(1) = 1, DIM(1) zeros are added at the bottom;
%    if LOC(2) = -1, DIM(2) zeros are added at the left;
%    if LOC(2) = 0, DIM(2)/2 zeros are added at the left and right;
%    if LOC(2) = 1, DIM(2) zeros are added at the right.
%  The default is LOC=[0,0].


if (nargin < 3), loc=[0,0]; end
s = size(a);
b = zeros(dim,class(a));

if (dim(1) == 0), dim(1)=s(1); end
if (dim(2) == 0), dim(2)=s(2); end

if (dim(1) < s(1))
  error('Size of output dim 1 (=%d) must be >= than input (=%d)', dim(1), s(1));
end
if (dim(2) < s(2))
  error('Size of output dim 2 (=%d) must be >= than input (=%d)', dim(2), s(2));
end

m = 1;  % pad bottom
k = 1;  % pad right
if (loc(1) == -1), m=dim(1)-s(1)+1; end  % pad top
if (loc(2) == -1), k=dim(2)-s(2)+1; end  % pad left
if (loc(1) == 0), m=ceil((dim(1)-s(1))/2+1); end  % pad symmetric
if (loc(2) == 0), k=ceil((dim(2)-s(2))/2+1); end  % pad symmetric

n = m + s(1) - 1;
l = k + s(2) - 1;

b(m:n, k:l) = a;
