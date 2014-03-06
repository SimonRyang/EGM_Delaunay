subroutine delaunay_int (xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash)

use comtypes
use comparam
use comvar
use comfunc
use delaunay_func

implicit none

real(prec),intent(out)::cdash,kdash,vxdash,vhdash,vfdash
real(prec)::xdash,hdash

real(prec)::w1,w2,w3,x1,x2,x3,h1,h2,h3

itri = save_tri

call walkt2 ( xdash, hdash, ntri, grid, til, tnbr, itri, iedg, ierr )

if (iedg .lt. 0) then
    print*, 'not on convex hull'
    convex = .false. 
end if

if (ierr .gt. 0) then
    print*, 'point not in triangle set- ', xc,hc,it
    ierr = 0
    call triangulation
    call walkt2 ( xdash, hdash, ntri, grid, til, tnbr, itri, iedg, ierr )
    if (ierr .gt. 0) then
        print*, 'point not in triangle set- again!!!! ', xc,hc,it
    end if
    save_tri = itri
end if

save_tri = itri

x1 = grid(1,til(1,itri))
x2 = grid(1,til(2,itri))
x3 = grid(1,til(3,itri))
h1 = grid(2,til(1,itri))
h2 = grid(2,til(2,itri))
h3 = grid(2,til(3,itri))

w2 = ((h3-h1)*(xdash-x3)+(x1-x3)*(hdash-h3))/((h2-h3)*(x1-x3)+(x3-x2)*(h1-h3))
w1 = ((h2-h3)*(xdash-x3)+(x3-x2)*(hdash-h3))/((h2-h3)*(x1-x3)+(x3-x2)*(h1-h3))
w3 = 1-w1-w2

if (w1 .lt. 0 .or. w2 .lt. 0 .or. w3 .lt. 0) then
   print*, 'point not in triangle', xc,hc,it
end if

cdash = w1*cf(til(1,itri))+w2*cf(til(2,itri))+w3*cf(til(3,itri))
vxdash = func_MUc(cdash)

kdash = w1*kf(til(1,itri))+w2*kf(til(2,itri))+w3*kf(til(3,itri))
vhdash = (netw+(1/gamma)*kdash**alpha)*func_MUc(cdash)

vfdash = w1*vf(til(1,itri))+w2*vf(til(2,itri))+w3*vf(til(3,itri))

end subroutine delaunay_int