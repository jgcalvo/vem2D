function L2err = getL2Error(out,sol)
% getL2Error computes the L2 error for the numerical solution of Poisson
% problem
%
% SYNOPSIS: L2err = getL2Error(out,sol)
%
% INPUT: out: structure generated by vem2d.m with the numerical solution
%	     sol: function handle with the exact solution
%
% OUTPUT: L2err: numerical approximation for the L2 error
%

% AUTHOR: Juan G. Calvo and collaborators, 2021

proj  = out.proj;  u     = out.u;
mesh  = out.mesh;  polys = out.polys;

for E = 1:size(mesh.elems,1)
    verts    = mesh.verts(mesh.elems{E},:);
    [h,~,xb] = geomElement(verts);
    tri  = delaunay(verts(:,1),verts(:,2));
    coef = proj{E}*u(mesh.dof{E});
    val  = 0;
    for k = 1:size(tri,1)
        xx = verts(tri(k,:),:);
        B11 = xx(2,1)-xx(1,1);
        B12 = xx(3,1)-xx(1,1);
        B21 = xx(2,2)-xx(1,2);
        B22 = xx(3,2)-xx(1,2);
        detB   = B11*B22-B12*B21;
        x = @(xhat,yhat) B11*xhat+B12*yhat+xx(1,1);
        y = @(xhat,yhat) B21*xhat+B22*yhat+xx(1,2);
        fun = @(xhat,yhat) (evalProj(x(xhat,yhat),y(xhat,yhat),coef,polys,xb,h)-sol(x(xhat,yhat),y(xhat,yhat))).^2;
        val = val + detB*integral2(fun,0,1,0,@(x) 1-x);
    end
    L2err = sqrt(val);
end
end

function polyn = evalProj(x,y,coef,polys,xb,h)
polyn = 0;
for poly_id = 1:size(polys,1)
    polyn = polyn + coef(poly_id)*((x-xb(1))/h).^polys(poly_id,1).*((y-xb(2))/h).^polys(poly_id,2);
end
end