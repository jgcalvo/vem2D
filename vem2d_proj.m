function proj = vem2d_proj(verts,k)
% vem2D_proj is a simplification of vem2d that computes only the projection
% operator for a given element; see vem2d.m for further details.
%
% SYNOPSIS: proj = vem2d_proj(verts,k)
%
% INPUT: verts:	vertices of the element
%	        k:	degree for local spaces (k>1)
%
% OUTPUT: proj: projector for the element
%

% AUTHOR: Juan G. Calvo and collaborators, 2021

% preliminaries
quad  = quadInfoGL(k);
polys = createPolynomials(k);
% constants and arrays for data
nk2     = (k-1)*k/2;
nPolys  = (k+1)*(k+2)/2;
nVerts = size(verts,1);
[h,area,centroid]  = geomElement(verts);
[nodesBd,vectorBd] = bndryNodesVectors(verts,quad.nodes);
intMon = precomputeIntegrals(verts,centroid,h,area,k);
base     = (nodesBd-centroid)/h;
n_dof    = k*nVerts + nk2;
% Matrix D
D = ones(n_dof, nPolys);      % first column is equal to one
for poly_id = 2:nPolys        % nodal dof
    D(1:k*nVerts,poly_id)  = prod(base.^polys(poly_id,:),2);
end
for dof_id = k*nVerts+1:n_dof % moments
    poly = polys(dof_id-k*nVerts,:)+polys;
    D(dof_id,:) = intMon(sub2ind(size(intMon), poly(:,1)+1, poly(:,2)+1))/area;
end
% Matrix B
B = zeros(nPolys, n_dof); B(1,k*nVerts+1) = 1;
for poly_id = 2:nPolys        % non-constant polynomials
    poly_degree = polys(poly_id,:);
    g1 = prod(base.^[max(poly_degree(1)-1,0) poly_degree(2)],2);
    g2 = prod(base.^[poly_degree(1) max(poly_degree(2)-1,0)],2);
    gradM_a = poly_degree.*[g1 g2]/h;
    B(poly_id, 1:k*nVerts) = sum(gradM_a.*vectorBd,2).*[repmat(quad.weights(1),nVerts,1); repmat(quad.weights(2:end-1),nVerts,1)];
end
for dof_id = k*nVerts+1:n_dof
    varphi_j = polys(dof_id-k*nVerts,:);
    poly_degree = polys(2:end,:);
    B(2:end,dof_id) = - area/h^2 * ...
        (poly_degree(:,1).*(poly_degree(:,1)-1).*prod(poly_degree-[2,0] == varphi_j,2)+poly_degree(:,2).*(poly_degree(:,2)-1).*prod(poly_degree-[0,2] == varphi_j,2));
end
proj = (B*D)\B;
end