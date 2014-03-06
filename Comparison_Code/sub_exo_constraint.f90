subroutine sub_exo_constraint(xb)

use comtypes
use comvar
use comfunc
use comfuncform
use sing_dim

implicit none

! general
real(prec)::hdash,xdash,vxdash,vhdash,vfdash

! brent
real(prec)::fsige 
real(prec)::xa,xb
logical::succes

xa = log(epsi)
xb = xa+log(ddx)

call  zbrac(f_exocon,xa,xb,succes)
xb = zbrent(f_exocon,xa,xb,tolbre)

fsige = f_exocon(xb)
if (abs(fsige) .gt. errabs) then
   print*, 'EXGM: failure of zbrent - borrowing constraint', xc,hc,it
   pause
endif

xb=1/(1+exp(-xb))*(gridxex(xc)+netw*gridhex(hc))

contains

!---------------------------------------------------------------------------------------
function f_exocon(x) 

real(prec),intent(in)::x
real(prec)::f_exocon
real(prec)::inv,consump,cdash,kdash

inv = 1/(1+exp(-x))*(gridxex(xc)+netw*gridhex(hc))

hdash = (1-delta)*(gridhex(hc)+func_F(inv))
xdash = 0.0

! Interpolate derivative of value function
call sub_exopolation(xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash)

consump = gridxex(xc)+netw*gridhex(hc)-inv
f_exocon=func_MUc(consump) - betta*(1-delta)*func_MFk(inv)*(func_margsurv(hdash)*vfdash+func_surv(hdash)*vhdash)

end function f_exocon
!---------------------------------------------------------------------------------------

end subroutine sub_exo_constraint