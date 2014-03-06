subroutine sub_hy_euler(x)

use comtypes
use comparam
use comvar
use comfunc
use comfuncform
!use sing_dim
use multi_dim

implicit none

real(prec),intent(out)::x(2)
real(prec)::hdash,xdash,cdash,kdash,vxdash,vhdash,vfdash,ctrans

!! zbrent
!real(prec)::xa,xb,fval1
!logical::succes

integer:: n                                   ! dimension of equation system
logical intlj,reevalj,check                   ! maximum step in line search algorithm (choose<1)
real(prec),dimension(1,1)::df,fval2,xf

!----------------------------------------------------------------------
! Rootfinding with Brent's method (Bisection/Secant)
! Initial Guess
!xa=epsi
!xb=xa+ddx
!call zbrac(f_hybrent,xa,xb,succes)
!xb=zbrent(f_hybrent,xa,xb,tolbr)
!x(2)=xb
!fval1=f_hybrent(xb)
!if (abs(fval1)>errabs) then
!    print*, 'ERROR: One dimensional rootfinding does not work', xc,hc,it
!endif
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Rootfinding with Broyden's method (Quasi Newton Method)
! Setting for Broyden
intlj = .true.
reevalj = .true.
n = 1  

! Initial guess for Investment in Human Capital
if (it .lt. 10) then
   if (xc .eq. 1 .and. hc .eq. 1) then
      xf = savekfunhy(xc,hc)
   elseif (xc .eq. 1 .and. hc .ne. 1) then
      xf = savekfunhy(xc,hc-1)
   else   
      xf = savekfunhy(xc-1,hc)
   endif
else
   xf = log(kfunhy(xc,hc))
endif

! Solving the 1-dimensional Equation System           
call broyden(f_hybroyden,fval2,xf,n,df,intlj,reevalj,check,maxstp,tolbro)

x(2) = exp(xf(1,1))

xdash = (1+ret)*grida(xc)
hdash = (1-delta)*(gridhhy(hc)+func_F(x(2)))

! Checking if State Space is sufficiently large
if (hdash .gt. gridhhy(nh)) then 
   print*,'h out of bounds, increase health on hand grid', xc,hc,it
end if
if (hdash .lt. gridhhy(1)) then
   print*, 'hdash.lt.gridhhy(1)) how can this be', xc,hc,it
endif

! Interpolate derivative of value function
call sub_hypolation(xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash)

ctrans = betta*(1+ret)*func_surv(hdash)*vxdash
x(1) = func_invuc(ctrans)

savekfunhy(xc,hc) = log(x(2))


contains
!!---------------------------------------------------------------------------------------
!function f_hybrent(inv)
!
!real(prec),intent(in)::inv
!real(prec):: f_hybrent
!
!hdash = (1-delta)*(gridhhy(hc)+func_F(inv))
!xdash = (1+ret)*grida(xc)
!
!! Interpolate derivative of value function
!call sub_hypolation(xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash)
!
!f_hybrent=func_MFk(inv)-((1+ret)/(1-delta))*(func_surv(hdash)*vxdash/(func_margsurv(hdash)*vfdash+func_surv(hdash)*vhdash))
!end function f_hybrent
!!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
subroutine f_hybroyden(fval,coeffs,nnm)

integer,intent(in)::nnm
real(prec),intent(in)::coeffs(nnm)
real(prec),intent(out)::fval(nnm)

real(prec):: inv(nnm)

inv = exp(coeffs)

hdash = (1-delta)*(gridhhy(hc)+func_F(inv(1)))
xdash = (1+ret)*grida(xc)

! Interpolate derivative of value function
call sub_hypolation(xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash)

fval(1)=func_MFk(inv(1))-((1+ret)/(1-delta))*(func_surv(hdash)*vxdash/(func_margsurv(hdash)*vfdash+func_surv(hdash)*vhdash))

end subroutine f_hybroyden
!---------------------------------------------------------------------------------------

end subroutine sub_hy_euler