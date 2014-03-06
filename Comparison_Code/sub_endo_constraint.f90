subroutine sub_endo_constraint(xb)

use comtypes
use comvar
use comfunc
use comfuncform
use sing_dim

implicit none

! general
real(prec)::hdash,xdash,cdash,kdash,vxdash,vhdash,vfdash

! brent
real(prec)::fsige
real(prec)::xa,xb
logical::succes

hdash = (1-delta)*gridz(hc)      ! human capital tomorrow
xdash = epsi                     ! assets tomorrow

! Interpolate derivative of value function
call delaunay_int(xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash)

xa = log(epsi)
xb = xa+log(ddx)

call  zbrac(f_endocon,xa,xb,succes)
xb = zbrent(f_endocon,xa,xb,tolbre)

fsige = f_endocon(xb)
if (abs(fsige) .gt. errabs) then
   print*, 'ENDGM: failure of zbrent - borrowing constraint', xc,hc,it
   pause
end if

xb = exp(xb)

contains

!---------------------------------------------------------------------------------------
function f_endocon(x) 

real(prec),intent(in)::x
real(prec)::f_endocon
real(prec):: inv,consump

inv = exp(x)

consump = aux_grid(xc)+netw*(hdash/(1-delta)-func_F(inv))-inv

f_endocon = func_MUc(consump)-betta*(1-delta)*func_MFk(inv)*(func_margsurv(hdash)*vfdash+func_surv(hdash)*vhdash)

end function f_endocon
!---------------------------------------------------------------------------------------

end subroutine sub_endo_constraint