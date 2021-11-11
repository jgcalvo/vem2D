function [out,mesh] = vem2d(mesh,k,f)
% VEM2D computes the virtual element solution of a Poisson problem on a
%	 polygonal mesh with degree k and right hand side f
%
% SYNOPSIS: out = vem2d(mesh,k,f)
%
% INPUT: mesh:	mesh structure (nodes, elements and boundary nodes)
%	        f:	handle function for the rhs of Poisson equation
%	        k:	degree for local spaces (k>1)
%
% OUTPUT: out:  struct with the following fields:
%               u: numerical approximation for the solution
%               h: max diameter of mesh
%               mesh: mesh with extra information
%               proj: projector for each element
%               polys: list of monomials
%               A: stiffness matrix
%               M: mass matrix
%

% AUTHOR: Juan G. Calvo and collaborators, 2021

% preliminaries
mesh  = meshSetup(mesh,k);     % include data on mesh
quad  = quadInfoGL(k);         % GL quad info
polys = createPolynomials(k);  % exponents for monomial basis

% constants and arrays for data
nk2       = (k-1)*k/2;
nPolys    = (k+1)*(k+2)/2;
projector = cell(size(mesh.elems,1),1); % projections for error estimates
ndofEl    = k*cellfun(@numel,mesh.elems)+nk2;    % dof per element
sizeMat   = sum(ndofEl.^2);
sizeRHS   = sum(ndofEl);
aa  = nan(sizeMat,1); ii   = nan(sizeMat,1);     % allocate space to
mm  = nan(sizeMat,1); jj   = nan(sizeMat,1);     % build mass and stiffness
rhs = nan(sizeRHS,1); irhs = nan(sizeRHS,1);     % matrices, and rhs
ptMat = 0;            ptRHS = 0;                 % pointers

% maximum of diameters
hmax = 0;

% iterate over elements
for E = 1:size(mesh.elems,1)
    % read vertices for element
    verts  = mesh.verts(mesh.elems{E},:);
    nVerts = size(verts,1);
    % dimension of local space
    n_dof    = k*nVerts + nk2;
    % compute diameter, area and centroid
    [h,area,centroid]  = geomElement(verts);
    % update maximum for diameter
    hmax = max(hmax,h);
    % precompute integral of monomials and vector on boundary
    intMon             = precomputeIntegrals(verts,centroid,h,area,k);
    [nodesBd,vectorBd] = bndryNodesVectors(verts,quad.nodes);
    % base of monomials
    base     = (nodesBd-centroid)/h;
    
    % assemble matrix D
    D = ones(n_dof, nPolys);      % first column is one
    for poly_id = 2:nPolys        % nodal dof
        D(1:k*nVerts,poly_id) = prod(base.^polys(poly_id,:),2);
    end
    for dof_id = k*nVerts+1:n_dof % moments
        poly = polys(dof_id-k*nVerts,:)+polys;
        D(dof_id,:) = intMon(sub2ind(size(intMon), poly(:,1)+1, poly(:,2)+1))/area;
    end
    
    % assemble matrix B
    B = zeros(nPolys, n_dof); B(1,k*nVerts+1) = 1;
    for poly_id = 2:nPolys              % loop for non-constant polynomials
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
    
    % local stiffness matrix
    G = B*D;
    PiStar = mypinv(G,B); %G\B;
    G(1, :) = 0;
    I = eye(n_dof);
    A = PiStar'*G*PiStar + (I-D*PiStar)'*(I-D*PiStar);
    projector{E} = PiStar;
    
    % mass matrix
    H = zeros(nPolys,nPolys);            % matrix H
    for i = 1:nPolys
        poly  = polys(i,:)+polys;
        H(i,:)= intMon(sub2ind(size(intMon),poly(:,1)+1,poly(:,2)+1));
    end
    C = zeros(nPolys,n_dof);             % matrix C
    C(1:nk2,k*nVerts+(1:nk2)) = area*eye(nk2);
    Ctem = H*PiStar;
    C(nk2+1:nPolys,:) = Ctem(nk2+1:nPolys,:);
    PiStar0 = mypinv(H,C);               % L2 projections Pi = H\C
    Pi0 = D*PiStar0;                     % H is ill-conditioned 
    M = C'*(PiStar0)+area*(I-Pi0)'*(I-Pi0);  % mass matrix
    
    % rhs
    int_fm = zeros(nPolys,1);            % compute (f,Pi*m_alpha)
    tri = delaunay(verts(:,1),verts(:,2));
    for j = 1:nPolys
        fun = @(x,y) (((x-centroid(1))./h).^polys(j,1)).*(((y-centroid(2))./h).^polys(j,2)).*f(x,y);
        int_fm(j) = integrate(fun,tri,verts);
    end
    rhsLoc = PiStar0'*int_fm;
    
    % global ordering - store local info
    [dofi,dofj] = meshgrid(mesh.dof{E},mesh.dof{E});
    ii(ptMat+1:ptMat+n_dof^2) = dofi(:);     % global index
    jj(ptMat+1:ptMat+n_dof^2) = dofj(:);     % global index
    mm(ptMat+1:ptMat+n_dof^2) = M(:);        % mass matrix
    aa(ptMat+1:ptMat+n_dof^2) = A(:);        % stiffness matrix
    rhs(ptRHS+1:ptRHS+n_dof)  = rhsLoc;      % rhs
    irhs(ptRHS+1:ptRHS+n_dof) = mesh.dof{E}; % index for rhs
    ptRHS = ptRHS + n_dof;
    ptMat = ptMat + n_dof^2;
end

% assemble matrices and rhs
A = sparse(jj,ii,aa);
M = sparse(jj,ii,mm);
b = accumarray(irhs,rhs(:),[size(A,1) 1]);

% save struct with info
out = struct('h',hmax,'mesh',mesh,'proj',{projector},...
    'polys',polys,'M',M,'A',A,'b',b);
end

function PiStar0 = mypinv(H,C)
    [U,S,V] = svd(H);
    ind = (diag(S)<1e-12*S(1,1));       % remove singular values close to 0
    Stem = diag(S); Stem(ind) = inf;
    PiStar0 = V*(diag(1./Stem)*(U'*C));
end