README file to describe the set of FORTRAN 90 programs that accompany the paper
'Endogenous Grids in Higher Dimensions: Delaunay Interpolation and Hybrid Methods'
by
Alexander Ludwig, ludwig@wiso.uni-koeln.de and
Matthias Schön,  m.schoen@wiso.uni-koeln.de

Copyright 2013 Alexander Ludwig and Matthias Schön
============================================================================

This replication package is available upon request.

The codes may be used and reproduced for educational and research purposes
freely, if all of the following three conditions are respected:
1. the above-mentioned paper is properly cited in all work and publications
   that rely at least partially on the methods implemented in these codes.
2. the copyright notice is reproduced with them,
3. the codes are not included in software sold for profit.
If your intended use of these codes stems from commercial, governmental 
or institutional interest, please contact one of the authors.

Make sure to extract the main program and all other programs contained
in this replication .zip archive into the same directory.

In the following we briefly describe the purpose of each of the subprograms
and provide information that might be helpful to those users who want to
apply these programs for solving problems of similar type.

Note that the programs contain plenty of comments and that the variable
names in the code match closely the notation in the paper. An exception is 
assets which are denoted by a in the paper and x in the code. The paper itself 
is to be considered the main source of documentation for understanding the method.

Settings:
Modifying algorithm and specifying parameters is done in comparam.f90
As the main function includes all three methods described in the paper one 
may only use a single one by setting the parameter ALGORITHM to the 
corresponding value. By choosing INFINITE one can either use the finite or 
the infinite horizon setting. There are other settings further explained by
comments in the code and all parameters are set according to the paper. The
exact grid size can also be varied in comparam.f90. 

Start
Running the main program comparison_main.f90 will execute the complete code
for given settings.

EXGM 
Main file is sub_exo.f90 and sub_exo_euler.f90. The EXGM interpolation is
in comfunc.f90. 

ENDGM 
Main file is sub_endo.f90 and sub_endo_euler.f90. The Delaunay
triangulation is done by triangulation.f90 the visibility walk is done by
walkt2 in delaunay_func.f90 and the Delaunay interpolation is done by
delauany_int.f90.

HEGM 
Main file is sub_hy.f90 and sub_hy_euler.f90. The HEGM interpolation is in
comfunc.f90.

Simulation
Using policy functions to simulate profiles to compute Euler Error and 
evaluate computational accuracy 

Output
Output is produced by sub_out.f90 which provides txt files for computational
speed, number of iteration till convergence and Euler equation errors.