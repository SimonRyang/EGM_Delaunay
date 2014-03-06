subroutine triangulation   

use comtypes
use comvar
use comparam
use delaunay_func

implicit none

integer f,k,j
   
f = 0
do j = 1,nh
   do k = 1,nx
      f = f+1
      grid(1,f) = gridxen(k,j)
      grid(2,f) = gridhen(k,j)
      cf(f) = cfunen(k,j)
      kf(f) = kfunen(k,j)
      vf(f) = vfunen(k,j)
      ind(f)= f
   end do   
end do

call dtris2 ( npt, npt, grid, ind, ntri, til, tnbr, stack, ierr )
if (ierr .gt. 0) then
   print*, 'failure in triangulation ', it    
end if

i_basis=.false.

end subroutine triangulation