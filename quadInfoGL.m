function quadrature = quadInfoGL(k)
% quadInfoGL computes the nodes and weights for the Gauss-Lobatto in [0,1]
%
% SYNOPSIS: quadrature = quadInfoGL(k)
%
% INPUT:  k:	degree
%
% OUTPUT: quadrature: struct with the fields 'nodes' and 'weights'
%

j = 1:k-2;
v = sqrt(j.*(j+2)./((2*j+1).*(2*j+3)));
A = diag(v,-1)+diag(v,1);
[V,D] = eig(A);
x  = diag(D);                           % quad nodes in [-1,1]
w  = V(1,:).^2*4/3;                     % quad weights
w1 = 1-sum(w./(1-x'))/2;                % weight initial node
w2 = 1-sum(w./(1+x'))/2;                % weight final node
quadrature.nodes   = ([-1; x; 1]+1)/2;  % translate to [0,1]
quadrature.weights = [w1, w./((1-x').*(1+x')), w2]'/2;
end