function val = integrate(f,tri,node)
% integrate computes the integral of a given function over a triangulation
% of a polygon.
%
% SYNOPSIS: val = integrate(f,tri,node)
%
% INPUT: f:	    handle function representing the integrand
%	     tri:	triangle connectivity 
%	     node:	coordinates of the vertices
%
% OUTPUT: val:  integral of f over the polygon
%

% AUTHOR: Juan G. Calvo and collaborators, 2021

val = 0;
for k = 1:size(tri,1)
    xx = node(tri(k,:),:);
    B11 = xx(2,1)-xx(1,1);
    B12 = xx(3,1)-xx(1,1);
    B21 = xx(2,2)-xx(1,2);
    B22 = xx(3,2)-xx(1,2);
    detB = B11*B22-B12*B21;
    x = @(xhat,yhat) B11*xhat+B12*yhat+xx(1,1);
    y = @(xhat,yhat) B21*xhat+B22*yhat+xx(1,2);
    h = @(xhat,yhat) f(x(xhat,yhat),y(xhat,yhat))*detB;
	val = val + integral2(h,0,1,0,@(x) 1-x);
end
end