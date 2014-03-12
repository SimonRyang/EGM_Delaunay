module comparam 

use comtypes

implicit none

! Problem Presettings     
integer,parameter::algorithm= 4               ! Execute: 1=EXGM; 2=ENDGM; 3=HYBRID; 4=All methods

logical,parameter::infinite = .true.         ! True: Infinite horizon setting; False: Finite horizon setting ending after nit iterations
integer,parameter::nit      = 999             ! number of maximal iteration steps/ number of years in finite setting
integer,parameter::age      = 50             ! number of maximal iteration steps/ number of years in finite setting

logical,parameter::smart    = .true.         ! 1: triangulation only in the first 50 iteration or if needed; 0:triangulation in every iteration

real(prec),parameter::epsi  = 1.0e-10         ! small number

! I Parameters
! Prices
real(prec),parameter::ret   = 0.03d0          ! Interest rate on assets
real(prec),parameter::netw  = 0.1d0           ! Net wage rate

! Preferences
real(prec),parameter::rhho  = 0.04d0          ! Time preference
real(prec),parameter::betta = 1.0/(1.0+rhho)  ! Annual discount factor
real(prec),parameter::tetta = 0.5d0		      ! Coefficient of relative risk aversion

! Human Capital Production
real(prec),parameter::alpha = 0.65d0          ! Elasticity of health production function !0.65
real(prec),parameter::gamma = 1.0d0           ! Production function
real(prec),parameter::delta = 0.05d0          ! Health depreciation

! Survival Probability
real(prec),parameter::phi   = 0.5             ! Adjustment parameter for survival rate  

! II Grids
! Size of State Space
integer,parameter::nh  = 300                  ! Points in human capital grid (exogenous/endogenous)
integer,parameter::nx  = nh    		          ! Points in asset grid (exogenous/endogenous)
integer,parameter::adx = 2                    ! Points in auxiliar grid - saving zero
integer,parameter::na  = nx-adx  	          ! Points in saving grid (exogenous)
integer,parameter::nz  = nh                   ! Points in human capital grid (exogenous)
integer,parameter::npt = nx*nh                ! Total Number of Grid Points

! Length of State Space
real(prec),parameter::mina = 0.0              ! min. value of savings
real(prec),parameter::maxa = 500.0            ! max. value of savings
real(prec),parameter::minz = 1.0              ! min. value of human capital savings
real(prec),parameter::maxz = 100.0            ! max. value of human capital savings
real(prec),parameter::minh = 0.95             ! min. value of human capital
real(prec),parameter::maxh = 100.0            ! max. value of human capital
real(prec),parameter::minx = 0.0              ! min. value of assets
real(prec),parameter::maxx = 500.0            ! max. value of assets

! Curvature of Grids
real(prec),parameter::curv = 3.0              ! Curvature of grids
real(prec),parameter::curv2= 3.0              ! Curvature of auxiliar grid

real(prec),parameter::crit = 1.0D-6           ! stopping criterion on policy functions

! III Numerical Solver
! Zbrent/brac
real(prec),parameter::tolbre = 1.0e-12        ! tolerance level of zbrent
real(prec),parameter::ddx    = 1.0e-4         ! range of first bracket in zbrac
real(prec),parameter::errabs = 1.0e-6         ! absolute error level 

! Broyden
real(prec),parameter::tolbro = 1.0e-12        ! abolute deviation of function values
real(prec),parameter::maxstp = 2.0            ! 

!V Simulation
real(prec),parameter::x0min = 100.0           ! Lower bound for initial assets for simulation
real(prec),parameter::x0max = 400.0           ! Upper bound for initial assets
real(prec),parameter::h0min = 40.0            ! Lower bound for initial human capital
real(prec),parameter::h0max = 80.0            ! Upper Bound for initial human capital
integer,parameter   ::nsim  = 10              ! number of simulations
  
end module comparam