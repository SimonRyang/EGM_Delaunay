subroutine sub_hy_constraint(xb)

use comtypes
use comvar
use comfunc
use comfuncform
use sing_dim

implicit none

! General
real(prec)::xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash

! Brent
real(prec)::fsige
real(prec)::xa,xb
logical::succes

xdash = 0.0

xa = log(epsi)
xb = xa+log(ddx)

call  zbrac(f_hycon,xa,xb,succes)
xb = zbrent(f_hycon,xa,xb,tolbre)

fsige = f_hycon(xb)
if (abs(fsige) .gt. errabs) then
   print*, 'HEGM: failure of zbrent - borrowing constraint', xc,hc,it
   pause
endif

xb = 1/(1+exp(-xb))*(aux_grid(xc)+netw*gridhhy(hc))

contains

!---------------------------------------------------------------------------------------
function f_hycon(x)

real(prec),intent(in)::x
real(prec)::f_hycon
real(prec)::inv,consump

inv = 1/(1+exp(-x))*(aux_grid(xc)+netw*gridhhy(hc))

hdash = (1-delta)*(gridhhy(hc)+func_F(inv))

! Interpolate derivative of value function
call sub_hypolation(xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash)

consump = aux_grid(xc)+netw*gridhhy(hc)-inv

f_hycon = func_MUc(consump)-betta*(1-delta)*func_MFk(inv)*(func_margsurv(hdash)*vfdash+func_surv(hdash)*vhdash)

end function f_hycon
!---------------------------------------------------------------------------------------

end subroutine sub_hy_constraint