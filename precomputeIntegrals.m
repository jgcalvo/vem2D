function intMon = precomputeIntegrals(verts,centroid,h,area,k)
% precomputeIntegrals computes integrals over an element 
% of the scaled monomials ((x-centroid)/h)^alpha, where alpha is a
% multiindex of degree at most 2k by using divergence theorem
%
% SYNOPSIS: intMon = precomputeIntegrals(verts,centroid,h,area,k)
%
% INPUT: verts:	   coordinates of the vertices (ordered counterclockwise)
%                  of the vertices of the element
%	     centroid: coordinates of the centroid of the element
%	     h:        diameter of the element
%        area:     area of the element
%        k:        polynomial degree
%
% OUTPUT: intMon:  (2k+1) x (2k+1) matrix with the value of integrals.
%                  Entry intMon(alpha1+1,alpha2+1) corresponds to the
%                  exponent alpha = (alpha1,alpha2).
%

% AUTHOR: Juan G. Calvo and collaborators, 2021

x  = verts(:,1); x(end+1) = x(1);
y  = verts(:,2); y(end+1) = y(1);
xb = centroid(1);   yb = centroid(2);
quadrature   = quadInfoGL(k+1);        % GL for degree 2k
polys        = createPolynomials(2*k); % monomials of deg 2k
intMon = nan(2*k+1,2*k+1);             % preallocate matrix
intMon(1,1) = area;
for index_monomial = 2:(2*k+1)*(k+1)   % integral for monomials
   int = 0;
   poly = polys(index_monomial,:);
   for i = 1:length(x)-1
    fun = @(t) (((x(i).*(1-t)+x(i+1).*t-xb)/h).^(poly(1)+1))...
       .*(((y(i).*(1-t)+y(i+1).*t-yb)/h).^(poly(2))).*(y(i+1)-y(i));
    int = int + sum(fun(quadrature.nodes).*quadrature.weights);
   end
   intMon(poly(1)+1,poly(2)+1) = h/(poly(1)+1)*int;
end
end