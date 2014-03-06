subroutine sub_exo

use comtypes
use comparam
use comvar
use comfunc, Only: sub_exopolation
use comfuncform
use clock

implicit none

real(prec):: xb,hprime
real(prec),dimension(2):: x
real(prec),dimension(nx,nh):: diffc,diffk

!*********** Exogenous is human capital and assets *********************************
method = 1

! Pre-Settings
savecfunex(1,1) = epsi
savekfunex(1,1) = epsi

!------------------------------------------------------------------------------------------------------------------

! Policy functions in last period
do hc = 1,nh
    do xc = 1,nx
        cfunex(xc,hc) = gridxex(xc)+netw*gridhex(hc)
    end do
end do
kfunex = 0.0

! Compute Value of last period for given consumpiton
do hc = 1,nh
    do xc = 1,nx
        vfunex(xc,hc) = func_U(cfunex(xc,hc))
    end do
end do

!------------------------------------------------------------------------------------------------------------------

! Iteration over the FOC using Exogenous Method / Broyden Algorithm
print*, '*********** Exogenous Method ************* '
call tic ! Start Time measurement

do it = 1,nit

    calt = cfunex
    kalt = kfunex 
    valt = vfunex

    marker = .false.

    ! Computing policy functions   
    do hc = 1,nh
       do xc = 1,nx
          call sub_exo_euler(x)

          ! If the borrowing constraint is violated solve onedimensional equation system
          if (marker(xc,hc)) then          
              call sub_exo_constraint(xb)
              x(2) = xb
              x(1) = gridxex(xc)+netw*gridhex(hc)-x(2) 
              if (x(1) .lt. 0.0 .or. x(2) .lt. 0.0) then 
                  print*, 'c or k lower than zero', xc,hc,it 
                  pause
              end if                                                 
          endif

          cfunex(xc,hc) = x(1)
          kfunex(xc,hc) = x(2)          
       end do ! xc       
    end do ! hc
        
    do hc = 1,nh
        do xc = 1,nx
           ! Update Value Function
           hprime = (gridhex(hc)+func_F(kfunex(xc,hc)))*(1-delta)
           vfunex(xc,hc) = func_U(cfunex(xc,hc))+betta*func_surv(hprime)*vfunex(xc,hc)
   
           ! Compute difference between old and new policy function
           diffc(xc,hc) = calt(xc,hc)-cfunex(xc,hc)
           diffk(xc,hc) = kalt(xc,hc)-kfunex(xc,hc)                  
        end do ! xc
    end do !hc
    
    ! Compute Maximum Distance between Policy Functions   
    criterion_c_ex(it) = maxval(abs(diffc(:,:)))
    criterion_k_ex(it) = maxval(abs(diffk(:,:)))  
    if (.not. infinite .and. it.le.age) then
        csimex(:,:,it) = cfunex(:,:) 
        ksimex(:,:,it) = kfunex(:,:)
        vsimex(:,:,it) = vfunex(:,:) 
    end if
    
    if (infinite) then
        if (criterion_c_ex(it) .lt. crit .and. criterion_c_ex(it) .lt. crit) then
            print* , 'EXGM finished after iteration step ', it   
            saveit(method) = it 
            exit
        end if 
    end if
end do ! Iteration 

if (infinite) then
    csimex(:,:,1) = cfunex(:,:) 
    ksimex(:,:,1) = kfunex(:,:)
    vsimex(:,:,1) = vfunex(:,:) 
end if 

call toc ! Stop time measurement
if (it .gt. nit .and. infinite)  then
    print* , 'EXGM not converged'
    saveit(method) = it
end if
timeins(method) = time_method(3,method)*60+time_method(4,method)
 
end subroutine sub_exo