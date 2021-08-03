function [nodesBd,vectorBd] = bndryNodesVectors(verts,quadNode)
% bndryNodesVectors computes vertices and vectors on the boundary of an
% element
%
% SYNOPSIS: [nodesBd,vectorBd] = bndryNodesVectors(verts,quadNode)
%
% INPUT: verts:	   coordinates of the vertices of the element
%        quadNode: Gauss-Lobatto quadrature nodes
%
% OUTPUT: nodesBd:  coordinates of the nodes related to boundary degrees of
%                   freedom, ordered counterclockwise. Vertices of the
%                   polygon are included first and then internal nodes on
%                   edges
%		  vectorBd: normal and average vectors for each nodal dof
%

% AUTHOR: Juan G. Calvo and collaborators, 2021


k = numel(quadNode)-1;   % degree
nVerts = size(verts,1);  % number of vertices
% dof on boundary (vertices and GL nodes)
initPts = repelem(verts,k-1,1);
direc   = repelem(diff([verts; verts(1,:)]),numel(quadNode)-2,1);
nodesBd = [verts; initPts+repmat(quadNode(2:end-1),nVerts,1).*direc];
% required vectors on boundary
tangV  = diff([verts; verts(1,:)]);     % tangential vectors
normal = [tangV(:,2), -tangV(:,1)];     % normal vectors
avNorm = normal+normal([end,1:end-1],:);% average normal
% vectors: average on vertices, normal on edges
vectorBd = [avNorm;reshape(reshape(repmat(normal,1,k-1)',1,2*(k-1)*nVerts),2,[])'];
end
