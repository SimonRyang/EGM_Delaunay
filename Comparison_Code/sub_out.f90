subroutine sub_out

use comtypes
use comparam
use comvar
implicit none

integer::gc 
real:: gridpoint
! Convergence Criterion
open (unit = 1, file = 'criterion_c.txt', status = 'old',position='append')
do gc = 1,nit
   write (1,'(3f30.15)') criterion_c_ex(gc), criterion_c_en(gc), criterion_c_hy(gc)
end do
close (1)

open (unit = 2, file = 'criterion_k.txt', status = 'old',position='append')
do gc = 1,nit
   write (2,'(3f30.15)') criterion_k_ex(gc), criterion_k_en(nit+1-gc), criterion_k_hy(gc)
end do
close (2)

! Iteration until convegence
if (infinite) then
    open (unit = 3, file = 'iterations.txt', status = 'old',position='append')!
    write (3, '(3i)') saveit(1), saveit(2), saveit(3)
    close (3)
end if

! Time measure
open (unit = 6, file = 'time.txt', status = 'old')
do gc = 1,4
      write (6, '(3f30.15)') time_method(gc,1),time_method(gc,2),time_method(gc,3)
end do
close(6)

gridpoint=real(nh)
open (unit = 7, file = 'timeins.txt', status = 'old',position='append')!,position='append'
write (7, '(4f30.15)') gridpoint, timeins(1), timeins(2), timeins(3)
close (7)

!open (unit = 8, file = 'error_max.txt', status = 'old',position='append')!,position='append'
!write (8, '(6f30.15)') mcx, mkx, mcn, mkn, mcy, mky
!close (8)
!
!open (unit = 9, file = 'error_av.txt', status = 'old',position='append')!,position='append'
!write (9, '(6f30.15)') acx, akx, acn, akn, acy, aky
!close (9)


end subroutine sub_out