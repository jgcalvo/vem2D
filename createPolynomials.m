function polys = createPolynomials(k)
% polys computes the exponent of the monomial basis in the natural order
% [0 0], [1 0], [0 1], [2 0], ... [0 k].
%
% SYNOPSIS: polys = createPolynomials(k)
%
% INPUT:   k:	degree for local spaces (k>1)
%
% OUTPUT: polys: (k+1)*(k+2)/2 x 2 matrix with list of exponents
%

% AUTHOR: Juan G. Calvo and collaborators, 2021

polys = zeros((k+1)*(k+2)/2,2);
pointer = 0;
for j = 0:k
    polys(pointer+1:pointer+j+1,:) = [j:-1:0; 0:j]';
    pointer = pointer + j + 1;
end