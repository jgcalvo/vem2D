% A numerical implementation for the high-order 2D VEM in MATLAB
% Version 1.0 27-Jul-2021
%
% Files
% 
% The directory contains the following files
% 
% README.txt            : this file
% vem2d.m               : main function for the virtual element method for 
%                         Poisson problem with homogeneous boundary conditions
% vem2d_proj.m          : similar as vem2d, but only computes the vem
%                         projector
% bndryNodesVectors.m   : compute GL nodes on the boundary of an element and 
%                         required normal and average vectors on the boundary
% createPolynomials.m   : create list of exponents for the monomial basis
% geomElement.m         : compute diameter, area and centroid of an element
% getL2Error.m          : approximate the L2 error of the numerical solution
% integrate.m           : auxiliary function that integrates a function over
%                         a general polygon (given by the union of triangles)
% meshSetup.m           : compute edge list, local to global ordering of dof
%                         and dof on the boundary
% precomputeIntegrals.m : precompute integrals over an element of the monomial
%                         basis by using divergence theorem
% quadInfoGL.m          : compute Gauss-Lobatto quadrature weights and nodes 
%                         for the interval [0,1]
% demo.m                : example with different options for mesh type and
%                         degree approximation
% 
% The 'meshes' subdirectory also contains the following files
% 
% hexagons_xx.mat      : data structure for five hexagonal meshes
% triangles_xx.mat     : data structure for five triangular meshes
% voronoi_xx.mat       : data structure for five Voronoi meshes
% 
% See the paper "A numerical implementation for the high-order 2D VEM in 
% MATLAB" for further explanations on the data structure.