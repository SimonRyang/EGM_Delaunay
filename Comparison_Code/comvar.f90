module comvar

use comtypes 
use comparam

implicit none

! Grids
! EXGM
real(prec),dimension(nx)::   gridxex                ! exogenous  assets grid
real(prec),dimension(nh)::   gridhex                ! exogenous  human capital grid
! ENDGM
real(prec),dimension(nx,nh)::gridxen                ! endogenous assets grid hand
real(prec),dimension(nx,nh)::gridhen                ! endogenous human capital grid
real(prec),dimension(nh)::   gridz                  ! exogenous  human capital savings grid

! ENDGM/HEGM
real(prec),dimension(nx)::   grida                  ! exogenous  savings grid

! HEGM
real(prec),dimension(nx,nh)::gridxhy                ! endogenous assets grid
real(prec),dimension(nh)::   gridhhy                ! exogenous  human capital grid

real(prec),dimension(adx+1)::aux_grid               ! for policy functions on borrowing constraint

! Value function
real(prec),dimension(nx,nh)::vfunex,vfunhy,vfunen   ! Value function

! Policy functions
real(prec),dimension(nx,nh)::cfunex,kfunex          ! Exogenous Method
real(prec),dimension(nx,nh)::cfunhy,kfunhy          ! Hybrid Method
real(prec),dimension(nx,nh)::cfunen,kfunen          ! Endogenous Method

real(prec),dimension(nx,nh)::calt,kalt,valt         ! Policy functions of previous step to compute derivatives
real(prec),dimension(nx,nh)::gridxalt,gridhalt      ! State space of previous step to compute derivatives    


real(prec),dimension(nx,nh,age)::csimex,ksimex,vsimex,csimen,ksimen,vsimen,csimhy,ksimhy,vsimhy,xsimhy
! Simulation for Delaunay 
real(prec),dimension(2,npt,age):: sav_grid 
real(prec),dimension(npt,age)::   sav_cf, sav_kf, sav_vf
integer,dimension(age):: sav_ntri  
integer,dimension(npt,age):: sav_ind                 
integer,dimension(3,2*npt,age):: sav_til, sav_tnbr

! Iteration Stopping Criterion
real(prec),dimension(nit)::  criterion_c_ex,criterion_k_ex  ! Exogenous Method
real(prec),dimension(nit)::  criterion_c_hy,criterion_k_hy  ! Hybrid Method
real(prec),dimension(nit)::  criterion_c_en,criterion_k_en  ! Endogenous Method

! Time measure
integer::method
real(prec),dimension(4,3)::  time_method             ! Computing time for algorithm
real(prec),dimension(3)::    timeins                 ! Computing time in seconds
integer,dimension(3)::       saveit                  ! number of iteration until convergence 

! Auxiliary variables for functions
integer:: xc,hc,it                                  ! Grid point marker
logical,dimension(nx,nh):: marker                   ! marker if borrowing constraint is binding in EXGM
logical convex

! Delaunay
integer:: ntri,ierr                                 
integer,dimension(npt)::ind                         
integer,dimension(3,2*npt)::til,tnbr                
integer,dimension(npt)::stack                       
integer:: itri,iedg                                 

integer:: save_tri                                  ! Starting triangle of new interpolation is final triangle form previous interpolation
real(prec),dimension(2,npt)::grid                   ! Grid of of next period/last iteration
real(prec),dimension(npt)::cf,kf,vf                 ! Policy functions of next period/last iteration   

!Broyden    
real(prec),dimension(nx,nh)::savecfunex,savekfunex  ! First guess two dimension
real(prec),dimension(nx,nh)::savekfunhy             ! First guess one dimension

! Basis
real(prec),dimension(2,2,2*npt):: a_inv
real(prec),dimension(2,1,2*npt):: const
integer:: count_nichtdoppel, count_doppel

logical,dimension(2*npt):: i_basis

end module comvar