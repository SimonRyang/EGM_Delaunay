subroutine sub_endo_euler(x)

use comtypes
use comparam
use comvar
use comfunc
use comfuncform

implicit none

real(prec),intent(out)::x(2)
real(prec)::hdash,xdash,cdash,kdash,vxdash,vhdash,vfdash,rhsc,rhsk

hdash = (1-delta)*gridz(hc)+epsi                 ! health on hand tomorrow
xdash = (1+ret)*(grida(xc))+epsi                 ! cash on hand tomorrow

! Delaunay interpolation of irregular grid to get derivatives of value function for xdash and hdash

call delaunay_int(xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash) 
!call delaunay_int_basis(xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash) 

rhsc = betta*(1+ret)*func_surv(hdash)*vxdash
x(1) = func_invuc(rhsc)
rhsk = ((1+ret)/(1-delta))*(func_surv(hdash)*vxdash/(func_margsurv(hdash)*vfdash+func_surv(hdash)*vhdash))
x(2) = func_invfk(rhsk)

end subroutine sub_endo_euler

