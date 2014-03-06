subroutine sub_hy

use comtypes
use comparam
use comvar
use comfunc
use comfuncform
use clock

implicit none

real(prec)::maxxhy,hprime
real(prec)::x_max_sav_min,x_min_sav_min,income,xb
real(prec),dimension(2)::x
real(prec),dimension(nx,nh):: diffc,diffk
integer count

!*********** Endogenous is health on hand, exogenous is cash on hand *********************************
method = 3

! Pre-Setting
savekfunhy = 1.0

!------------------------------------------------------------------------------------------------------------------

! Define initial grid for assets
maxxhy = grida(nx)*(1.0+ret)+epsi 
do hc = 1,nh
   gridxhy(:,hc) = func_makegrid(minx,maxxhy,nx,curv)
end do

! Policy functions in last period
do hc = 1,nh
    do xc = 1,nx
        cfunhy(xc,hc) = gridxhy(xc,hc)+netw*gridhhy(hc) 
    end do
end do  
kfunhy = 0.0

! Compute Value of last period for given consumpiton
do hc = 1,nh
    do xc = 1,nx
        vfunhy(xc,hc) = func_U(cfunhy(xc,hc))
    end do
end do

!------------------------------------------------------------------------------------------------------------------

! Iteration over the FOC using Hybrid Method / Broyden/Brac-Brent 
print*, '*********** Hybrid Method **************'
call tic ! Start Time measurement

do it = 1,nit
   
   calt = cfunhy
   kalt = kfunhy
   valt = vfunhy
   gridxalt = gridxhy
   
   ! Computing the policy functions
   do hc = 1,nh 
      ! Endogenous part of cash on hand grid
      do xc = adx+1,nx  
         call sub_hy_euler(x)
         cfunhy(xc,hc) = x(1)
         kfunhy(xc,hc) = x(2)
      end do
    
      ! Exogenous part of cash on hand grid (saving is zero)
      income = netw*gridhhy(hc)
      x_max_sav_min = grida(adx+1)-income+cfunhy(adx+1,hc)+kfunhy(adx+1,hc) ! max. cash on hand when optimal saving is zero  
      x_min_sav_min = min(-income/2.0,0.0)
      aux_grid(1:adx+1) = func_makegrid( x_min_sav_min,x_max_sav_min,adx+1,curv2)
      do xc = 1,adx
         call sub_hy_constraint(xb)
         kfunhy(xc,hc) = xb
         cfunhy(xc,hc) = aux_grid(xc)+netw*gridhhy(hc)-kfunhy(xc,hc)
      end do
   end do

   do hc = 1,nh
      do xc = 1,nx
         ! Update endogenous grids
         gridxhy(xc,hc)=grida(xc)-netw*gridhhy(hc)+cfunhy(xc,hc)+kfunhy(xc,hc)
  
         ! Update Value Function
         hprime=(1-delta)*(gridhhy(hc)+func_F(kfunhy(xc,hc)))
         vfunhy(xc,hc)=func_U(cfunhy(xc,hc))+betta*func_surv(hprime)*vfunhy(xc,hc)
   
         ! Compute difference between old and new policy function
         diffc(xc,hc)=calt(xc,hc)-cfunhy(xc,hc)
         diffk(xc,hc)=kalt(xc,hc)-kfunhy(xc,hc)
      end do ! xc
   end do ! hc
    
    ! Compute Maximum Distance between Policy Functions
    criterion_c_hy(it) = maxval(abs(diffc(:,:)))
    criterion_k_hy(it) = maxval(abs(diffk(:,:)))   
    
    if (.not. infinite .and. it.le.age) then    
        csimhy(:,:,it) = cfunhy(:,:) 
        ksimhy(:,:,it) = kfunhy(:,:)  
        xsimhy(:,:,it) = gridxhy(:,:)
        vsimhy(:,:,it) = vfunhy(:,:)
    end if
        
    if (infinite) then
        if (criterion_c_hy(it) .lt. crit .and. criterion_c_hy(it) .lt. crit) then
            print* , 'HEGM finished after iteration step ', it  
            saveit(method) = it 
            exit
        end if
    end if
end do ! Iteration 

if (infinite) then
    csimhy(:,:,1) = cfunhy(:,:) 
    ksimhy(:,:,1) = kfunhy(:,:)  
    xsimhy(:,:,1) = gridxhy(:,:)
    vsimhy(:,:,1) = vfunhy(:,:)
end if 

call toc ! Stop time measurement
if (it .gt. nit .and. infinite)  then
    print* , 'HEGM not converged'
    saveit(method) = it
end if
timeins(method) = time_method(3,method)*60+time_method(4,method)

end subroutine sub_hy