program comparison_main

use comparam

implicit none

! Establishing Grids
call sub_grids

print*, '*************** Solution ************* '
! Running different methods
if (algorithm .eq. 1 .or. algorithm .eq. 4) call sub_exo            ! Method with two exogenous grids

if (algorithm .eq. 2 .or. algorithm .eq. 4) call sub_endo           ! Method with two endogenous grids

if (algorithm .eq. 3 .or. algorithm .eq. 4) call sub_hy             ! Method with one exogenous and one endogenous grid

! Producing output

print*, '*************** Simulation ************* '
call sub_sim                ! Running Simulations to Compute Euler Errors

call sub_out                ! Save Output in txt.file
 
!pause   
    
end program comparison_main