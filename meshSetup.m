function mesh = meshSetup(mesh,k)
% meshSetup computes required fields on the mesh for the virtual element
% method
%
% SYNOPSIS: mesh = meshSetup(mesh,k)
%
% INPUT: mesh:	mesh structure (nodes, elements and boundary nodes)
%	        k:	degree for local spaces (k>1)
%
% OUTPUT: mesh: struct with the following new fields:
%               edges: list of edge connectivity on the mesh
%               dof:   local to global mapping of dof
%               bdDOF: global dof on the boundary of the domain
%

% AUTHOR: Juan G. Calvo and collaborators, 2021

nNodes = size(mesh.verts,1);           % global number of nodes
nVerts = cellfun(@numel,mesh.elems);   % vertices per element
% create array with all edges [e1 e2] (smaller index first)
e1 = cell2mat(mesh.elems);
e2 = cell2mat(cellfun(@(y)y([2:end,1]),mesh.elems,'UniformOutput',0));
[edges,index] = sort([e1 e2],2);
signEdges = index(:,2)-index(:,1);     % reversed edges
[mesh.edges, ~, indexEdges] = unique(edges,'rows');
numEdges = size(mesh.edges,1);         % global number of edges
% find edges on boundary
uniqueEdges = hist(indexEdges,1:size(mesh.edges,1))';
edgesBd = uniqueEdges==1;
edgesBd = nNodes+(find(edgesBd)-1)*(k-1)+repmat(1:k-1,sum(edgesBd),1);
mesh.bdDOF = sort([mesh.bndry; edgesBd(:)]);
% create array with global dof for each element
index = [0 cumsum(nVerts)'];
for j = 1:size(mesh.elems,1)
    range = index(j)+1:index(j+1);
    increment = repmat(1:k-1,nVerts(j),1); 
    % reverse order on edges that changed orientation
    sign  = signEdges(range);
    increment(sign<0,:) = fliplr(increment(sign<0,:));
    % find global dof for local vertices, edges and moments
    dofVert=mesh.elems{j};
    dofEdgs=(nNodes+(indexEdges(range)-1)*(k-1)+increment)';
    dofElem=nNodes+numEdges*(k-1)+(j-1)*k*(k-1)/2+(1:k*(k-1)/2)';
    mesh.dof{j} = [dofVert; dofEdgs(:); dofElem];
end
end