  ***************************************************************************
  * All the software  contained in this library  is protected by copyright. *
  * Permission  to use, copy, modify, and  distribute this software for any *
  * purpose without fee is hereby granted, provided that this entire notice *
  * is included  in all copies  of any software which is or includes a copy *
  * or modification  of this software  and in all copies  of the supporting *
  * documentation for such software.                                        *
  ***************************************************************************
  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
  ***************************************************************************
  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
  * ABOVE STATEMENT.                                                        *
  ***************************************************************************

   AUTHORS:

       Juan G. Calvo
       Universidad de Costa Rica, Costa Rica
       E-mail: juan.calvo@ucr.ac.cr

       César Herrera
       Universidad de Costa Rica, Costa Rica
       E-mail: cesar.herreragarro@ucr.ac.cr

       Ricardo Corrales-Barquero
       Universidad de Costa Rica, Costa Rica
       E-mail: ricardo.corralesbarquero@ucr.ac.cr

       Jorge Arroyo-Esquivel
       University of California Davis, United States
       E-mail: jarroyoe@ucdavis.edu

   REFERENCE:

       A numerical implementation for the high-order 2D VEM in MATLAB. Submitted.

   SOFTWARE REVISION DATE:

       V1.0, July 2021

   SOFTWARE LANGUAGE:

       MATLAB R2020b and later


============================================================================
SOFTWARE
============================================================================
This software provides all the functions needed to compute approximate
solutions to Poisson's problem using the virtual element method for a given
polynomial degree on polygonal meshes.

============================================================================
PACKAGE
============================================================================

The directory contains the following files

README.txt            : this file
vem2d.m               : main function for the virtual element method for 
                        Poisson problem with homogeneous boundary conditions
vem2d_proj.m          : similar as vem2d, but only computes the vem
                        projector
bndryNodesVectors.m   : compute GL nodes on the boundary of an element and 
                        required normal and average vectors on the boundary
createPolynomials.m   : create list of exponents for the monomial basis
geomElement.m         : compute diameter, area and centroid of an element
getL2Error.m          : approximate the L2 error of the numerical solution
integrate.m           : auxiliary function that integrates a function over
                        a general polygon (given by the union of triangles)
meshSetup.m           : compute edge list, local to global ordering of dof
                        and dof on the boundary
precomputeIntegrals.m : precompute integrals over an element of the monomial
                        basis by using divergence theorem
quadInfoGL.m          : compute Gauss-Lobatto quadrature weights and nodes 
                        for the interval [0,1]
demo.m                : example with different options for mesh type and
                        degree approximation

The 'meshes' subdirectory also contains the following files

hexagons_xx.mat      : data structure for five hexagonal meshes
triangles_xx.mat     : data structure for five triangular meshes
voronoi_xx.mat       : data structure for five Voronoi meshes

See the paper "A numerical implementation for the high-order 2D VEM in 
MATLAB" for further explanations on the data structure.

============================================================================
HOW TO INSTALL
============================================================================

Unpack vem2d.zip. The directory vem2d will be created with all the needed 
files inside. Add the directory vem2d to the path in Matlab or change the 
current working directory to vem2d.

============================================================================
HOW TO USE
============================================================================

Running the method simply involves calling the 'vem2d' function with the
following three arguments:

    1. the mesh given by a struct 'mesh'
    2. the polynomial degree for the local spaces
    3. a function handle for the right hand side of the PDE

See demo.m for different examples.

=============================================================================
BUG FIXES / UPDATES
=============================================================================

The most recent version of this package can be found on Github, at the
address https://github.com/jgcalvo.