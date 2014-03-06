subroutine sub_exo_euler (x)

use comtypes
use comparam
use comvar
use comfunc
use comfuncform 
use multi_dim

implicit none

real(prec),intent(out)::x(2)
real(prec)::xdash,hdash,vxdash,vhdash,vfdash

integer:: n                                   ! dimension of equation system
logical intlj,reevalj,check                   ! maximum step in line search algorithm (choose<1)
real(prec)::fval(2)
real(prec)::df(2,2)

! Setting for Broyden
intlj   = .true.
reevalj = .true.
n = 2                                         

! Initial Guess for Consumption and Investment in Human Capital
if (it .lt. 10) then
   if (xc .eq. 1 .and. hc .eq. 1) then
      x(1) = savecfunex(xc,hc)
      x(2) = savekfunex(xc,hc)
   elseif (xc .eq. 1 .and. hc .ne. 1) then
      x(1) = savecfunex(xc,hc-1)
      x(2) = savekfunex(xc,hc-1)
   else   
      x(1) = savecfunex(xc-1,hc)
      x(2) = savekfunex(xc-1,hc)
   endif
else
   x(1) = log(cfunex(xc,hc))
   x(2) = log(kfunex(xc,hc))
endif

! Solving the 2-dimensional Equation System
call broyden(f_exo,fval,x,n,df,intlj,reevalj,check,maxstp,tolbro)

x = exp(x)

! Checking if State Space is sufficiently large
xdash = (1+ret)*(gridxex(xc)+netw*gridhex(hc)-x(1)-x(2))
hdash = (1-delta)*(gridhex(hc)+func_F(x(2))) 

!if (hdash .gt. gridhex(nh)) then
!   print*,'h out of bounds, increase max human capital', xc,hc,it
!end if
!if (xdash .gt. gridxex(nx)) then
!   print*,'x out of bounds, increase max assets', xc,hc,it
!end if
!if (hdash .lt. gridhex(1)) then
!   print*,'h out of bounds, increase min human capital', xc,hc,it
!end if

! Is borrowing constraint violated
if (xdash.lt.0.0) then
    marker(xc,hc)=.true.
end if 

savecfunex(xc,hc)=log(x(1))
savekfunex(xc,hc)=log(x(2))

contains

!---------------------------------------------------------------------------------------
subroutine f_exo(fval,coeffs,nnm)

implicit none

integer,intent(in)::nnm
real(prec),intent(out)::fval(nnm)
real(prec),intent(in)::coeffs(nnm)

real(prec)::x1(2),cdash,kdash

x1 = exp(coeffs)

! assets and human capital tomorrow
xdash = (1 + ret) * (gridxex(xc)+netw*gridhex(hc)-x1(1)-x1(2))
hdash = (1-delta) * (gridhex(hc)+func_F(x1(2)))

! Interpolate derivative of value function
call sub_exopolation(xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash)

! Equation System (FOC)
fval(1) = func_MUc(x1(1)) - betta*(1+ret)*func_surv(hdash)*vxdash 
fval(2) = func_MFk(x1(2)) - ((1+ret)/(1-delta))*(func_surv(hdash)*vxdash/(func_margsurv(hdash)*vfdash+func_surv(hdash)*vhdash))

end subroutine f_exo
!---------------------------------------------------------------------------------------

end subroutine sub_exo_euler