subroutine sub_grids

use comtypes 
use comparam
use comvar
use comfunc, only: func_makegrid
implicit none

! EXGM
gridxex(1:nx) = func_makegrid(minx,maxx,nx,curv)      ! Grid for assets in EXGM
gridhex(1:nh) = func_makegrid(minh,maxh,nh,curv)      ! Grid for human capital in EXGM

! ENDGM
gridz(1:nh) = func_makegrid(minz,maxz,nh,curv)        ! Grid for human capital savings 

! ENDGM/HEGM
grida(adx+1:nx) = func_makegrid(mina,maxa,na,curv)    ! Grid for savings
grida(1:adx) = mina                                   ! Binding part of borrowing constraint savings set to zero   

! HEGM
gridhhy(1:nh) = func_makegrid(minh,maxh,nh,curv)      ! Grid for human capital in HEGM

end subroutine sub_grids