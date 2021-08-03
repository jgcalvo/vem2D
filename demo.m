function demo(varargin)
% demo includes two examples for the virtual element solution of a Poisson
% problem on a polygonal mesh
%
% SYNOPSIS: demo(opt,meshType)
%
% INPUT:  opt:	    scalar with two options:
%                   1: convergence example as a function of h and k for a
%                   small set of values of h and k (due to running times)
%                   2: plot solution for a particular value of h and k
%                   3: approximate harmonic functions without VEM for the
%                   unit square (meshType is not required)
%                   4: approximate harmonic functions without VEM for a
%                   L-shaped domain (meshType is not required)
%         meshType: type of mesh; use 'tri','hex' or 'vor' to use
%                   provided meshes for options 1 and 2
%
% EXAMPLE:
%         demo(1,'tri'),

% AUTHOR: Juan G. Calvo and collaborators, 2021

if(nargin == 2)
    opt = varargin{1};
    meshType = varargin{2};
else
    opt = varargin{1};
end

switch opt
    case 1 % convergence as a function of h and k
        % load mesh files for triangles; can also use voronoi or hexagons
        files = dir(['./meshes/' meshType '*']);
        % polynomial degree
        k = 2:2:6;
        % matrices to store data
        errL2  = nan(3,numel(k));
        diam   = nan(3,numel(k));
        % rhs and exact solution
        f   = @(x,y) 2*pi*pi*sin(pi*x).*sin(pi*y);
        sol = @(x,y) sin(pi*x).*sin(pi*y);
        
        for indFile = 1:3 % mesh loop
            mesh = load(['./meshes/' files(indFile).name]);
            mesh = mesh.mesh;
            for indK = 1:numel(k)    % degree loop
                out = vem2d(mesh,k(indK),f); % VEM main function
                % compute L2 error
                errL2(indFile,indK)  = getL2Error(out,sol);
                % store diameter
                diam(indFile,indK)   = out.h;
            end
        end
        
        % plot L2 error as fn of
        figure(1)
        indMesh = 2;
        loglog(k,errL2(indMesh,:))
        xlabel('k')
        ylabel('L2 Error')
        
        % plot L2 error as fn of h
        figure(2)
        for indK = 1:numel(k)
            loglog(diam(:,indK),errL2(:,indK)), hold on
            textLegend{indK} = ['k=' num2str(k(indK))];
        end
        xlabel('h')
        ylabel('L2 Error')
        legend(textLegend,'Location','NorthWest')
    case 2 % plot solution for a particular choice of h and k
        file = dir(['./meshes/' meshType '*_02.mat']);
        mesh = load(['./meshes/' file(1).name]);
        mesh = mesh.mesh;
        k = 4;
        % rhs and exact solution
        f   = @(x,y) 2*pi*pi*sin(pi*x).*sin(pi*y);
        sol = @(x,y) sin(pi*x).*sin(pi*y);
        out = vem2d(mesh,k,f); % VEM main function
        % plot numerical solution
        % we only consider values on vertices for simplicity
        valVertex = out.u(1:size(mesh.verts,1));
        elems = convertCells(mesh.elems);
        figure(1)
        patch('Faces', elems, 'Vertices', mesh.verts,'CData', valVertex,'FaceColor','interp');
        colorbar
        axis equal
        axis off
        % plot numerical error
        figure(2)
        error = abs(valVertex-sol(mesh.verts(:,1),mesh.verts(:,2)));
        %cData = cData/norm(cData,inf);
        patch('Faces', elems, 'Vertices', mesh.verts,'CData', error,'FaceColor','interp');
        colorbar
        axis equal
        axis off
    case {3, 4} % mesh free solver for approximating harmonic functions
        % store error
        error = zeros();  pt = 1;
        % degrees
        k = 2:9;
        % domain and exact solution
        if(opt == 3) % unit square
            verts = [0 0; 1 0; 1 1; 0 1];
            uex = @(x) 1*(exp(2*x(:,1))+exp(-2*x(:,1))).*sin(2*x(:,2));
        else         % L shaped domain
            verts = [0 0; 1 0; 1 1; -1 1; -1 -1; 0 -1];
            uex = @exactu;
        end
        % iterate for each degree
        for max_degree = k
            nnk = max_degree*(max_degree-1)/2;
            quadrature  = quadInfoGL(max_degree);
            proj = vem2d_proj(verts,max_degree);
            polys = createPolynomials(max_degree);
            [diam,area,centroid]  = geomElement(verts);
            [nodesBd,~] = bndryNodesVectors(verts,quadrature.nodes);
            intMon = precomputeIntegrals(verts,centroid,diam,area,max_degree);
            nPolys = size(polys,1);
            % matrix with entries \int_E \nabla m_i \nabla m_j
            Q = zeros(nPolys-1,nPolys-1);
            for j = 2:nPolys
                alpha = polys(j,:);
                for m = 2:size(polys,1)
                    beta  = polys(m,:);
                    tem = 0;
                    if(alpha(1)+beta(1)>1)
                        tem = tem + alpha(1)*beta(1)/diam^2*intMon(alpha(1)+beta(1)-1,alpha(2)+beta(2)+1);
                    end
                    if(alpha(2)+beta(2)>1)
                        tem = tem + alpha(2)*beta(2)/diam^2*intMon(alpha(1)+beta(1)+1,alpha(2)+beta(2)-1);
                    end
                    Q(j-1,m-1) = tem;
                end
            end
            D1 = uex(nodesBd);  % nodal dof
            P = proj(2:end,:);  % remove constant monomial
            mat = P'*Q*P;
            M22 = mat(end-nnk+1:end,end-nnk+1:end);
            M21 = mat(end-nnk+1:end,1:end-nnk);
            D2 = -M22\M21*D1;   % interior dof
            dof = [D1; D2];
            coef = proj*dof;    % coef of solution (in monomial basis)
            % evaluate at point [.5 .5]
            x = [.5 .5]; 
            M = ((x(1)-centroid(1))/diam).^polys(:,1).*((x(1)-centroid(2))/diam).^polys(:,2);
            error(pt) = abs(M'*coef-uex(x));
            pt = pt + 1;
        end
        loglog(k,error)
        xlabel('k')
        ylabel('Error at one point')
end
end

function elemsMatrix = convertCells(elems) % convert cell elems to matrix
maxVerts = max(cellfun(@numel,elems));
numElems = size(elems,1);
elemsMatrix = NaN(numElems,maxVerts);
for j = 1:numElems
    elemsMatrix(j,1:numel(elems{j})) = elems{j};
end
end

function u = exactu(p) % solution for L shaped domain
r = sqrt(sum(p.^2,2));
theta = atan2(p(:,2),p(:,1));
theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
u = r.^(2/3).*sin(2*theta/3);
end