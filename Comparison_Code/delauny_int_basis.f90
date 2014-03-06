subroutine delaunay_int_basis (xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash)

use comtypes
use comparam
use comvar
use comfunc
use delaunay_func

implicit none

real(prec),intent(out)::cdash,kdash,vxdash,vhdash,vfdash
real(prec)::xdash,hdash

real(prec),dimension(2,1):: B, T

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

if (i_basis(itri)==.false.) then
    call sub_basis
    i_basis(itri)=.true.
    count_nichtdoppel=count_nichtdoppel+1  
else
    count_doppel=count_doppel+1
end if

B(1,1) = xdash
B(2,1) = hdash
T = matmul(A_inv(:,:,itri),B)-CONST(:,:,itri)

cdash = cf(til(1,itri)) + T(1,1)*(cf(til(2,itri))-cf(til(1,itri))) + T(2,1)*(cf(til(3,itri))-cf(til(2,itri)))
vxdash = func_MUc(cdash)

kdash = kf(til(1,itri)) + T(1,1)*(kf(til(2,itri))-kf(til(1,itri))) + T(2,1)*(kf(til(3,itri))-kf(til(2,itri)))
vhdash = (netw+(1/gamma)*kdash**alpha)*func_MUc(cdash)

vfdash = vf(til(1,itri)) + T(1,1)*(vf(til(2,itri))-vf(til(1,itri))) + T(2,1)*(vf(til(3,itri))-vf(til(2,itri)))

end subroutine delaunay_int_basis