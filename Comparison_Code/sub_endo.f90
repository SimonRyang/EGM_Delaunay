subroutine sub_endo

use comtypes
use comparam
use comvar
use comfunc
use comfuncform
use clock

implicit none

real(prec)::maxhen,maxxen
real(prec)::x_max_sav_min,x_min_sav_min,income
real(prec)::xb
real(prec),dimension(2)::x
real(prec),dimension(nx,nh):: diffc,diffk
integer count

integer:: f,j,k

!*********** Endogenous is health on hand and cash on hand *********************************
method = 2

sav_grid = 0.0
sav_cf = 0.0
sav_kf = 0.0
sav_vf = 0.0
sav_ntri = 0.0
sav_ind  = 0.0        
sav_til = 0.0
sav_tnbr = 0.0

!------------------------------------------------------------------------------------------------------------------

! Define initial grids for assets and human capital
maxxen = grida(nx)*(1.0+ret)+2*epsi
maxhen = gridz(nz)*(1-delta)+epsi
do hc = 1,nh
   do xc = 1,nx
      gridxen(:,hc) = func_makegrid(minx,maxxen,nx,curv)
      gridhen(xc,:) = func_makegrid((1-delta)*minz,maxhen,nh,curv2)!
   end do
end do

! Policy functions in last period
cfunen = gridxen+netw*gridhen
kfunen = 0.0

! Compute Value of last period for given consumption
do hc = 1,nh
    do xc = 1,nx
        vfunen(xc,hc) = func_U(cfunen(xc,hc))
    end do
end do


!------------------------------------------------------------------------------------------------------------------

! Iteration over the FOC using Endogenous Method / Delaunay
print*, '*********** Endogenous Method ************* '
call tic ! Start Time measurement

do it = 1,nit 
    gridxalt = gridxen 
    gridhalt = gridhen
    calt     = cfunen
    kalt     = kfunen
    valt     = vfunen

   ! Delaunay Triangulation of Irregular Grid
   i_basis=.false.
   if (smart) then    
      if (it .lt. 50) then ! .or. mod(it,100).eq.0
         call triangulation
         save_tri = 1 
      end if
      
      f = 0
      do j = 1,nh
         do k = 1,nx
            f = f+1
            grid(1,f) = gridxen(k,j)
            grid(2,f) = gridhen(k,j)
            
            cf(f) = cfunen(k,j)
            kf(f) = kfunen(k,j)
            vf(f) = vfunen(k,j)
         end do
      end do
   else
       call triangulation
       save_tri = 1 
   end if 
          
  ! Computing Policy functions
   do hc = 1,nh                                  
        ! Endogenous part of asset grid
        do xc = adx+1,nx
            call sub_endo_euler(x)
            cfunen(xc,hc) = x(1)
            kfunen(xc,hc) = x(2)
        end do

        ! Exogenous part of asset grid (saving is zero)
        income = netw*(gridz(hc)-func_F(kfunen(adx+1,hc)))
        x_max_sav_min = grida(adx+1)-income+cfunen(adx+1,hc)+kfunen(adx+1,hc) ! max. cash on hand when optimal saving is zero
        x_min_sav_min = min(-income/2.0,0.0)
        aux_grid(1:adx+1) = func_makegrid(x_min_sav_min,x_max_sav_min,adx+1,curv2)    ! auxiliar grid for computing policy function at binding borrowing constraint
        do xc = 1,adx
           call sub_endo_constraint(xb)
           kfunen(xc,hc) = xb
           cfunen(xc,hc) = aux_grid(xc)+netw*(gridz(hc)-func_F(kfunen(xc,hc)))-kfunen(xc,hc)
        end do
   end do 
   
   do hc = 1,nh
      do xc = 1,nx
         ! Update endogenous grids
         gridhen(xc,hc) = gridz(hc)-func_F(kfunen(xc,hc)) 
         gridxen(xc,hc) = grida(xc)-netw*gridhen(xc,hc)+cfunen(xc,hc)+kfunen(xc,hc)
         
         ! Update Value Function
         vfunen(xc,hc) = func_U(cfunen(xc,hc))+betta*func_surv((1-delta)*gridz(hc))*vfunen(xc,hc)
                  
         ! Compute difference between old and new policy function
         diffc(xc,hc) = calt(xc,hc)-cfunen(xc,hc)
         diffk(xc,hc) = kalt(xc,hc)-kfunen(xc,hc)
      end do
   end do
    
    ! Compute Maximum Distance between Policy Functions
    criterion_c_en(it) = maxval(abs(diffc(:,:)))
    criterion_k_en(it) = maxval(abs(diffk(:,:)))
    
    if (.not. infinite .and. it.le.age) then 
            sav_grid(1,:,it) = grid(1,:) 
            sav_grid(2,:,it) = grid(2,:) 
            sav_cf(:,it) = cf(:) 
            sav_kf(:,it) = kf(:) 
            sav_vf(:,it) = vf(:)             

            sav_ntri(it) = ntri                                
            sav_ind(:,it) = ind(:)                         
            sav_til(:,:,it) = til(:,:)                
            sav_tnbr(:,:,it) = tnbr(:,:)  
        end if  
        
    if (infinite) then
        if (criterion_c_en(it) .lt. crit .and. criterion_c_en(it) .lt. crit) then
            print* , 'ENDGM finished after iteration step ', it   
            saveit(method) = it 
            exit
        end if
    end if   
end do ! Iteration 

call toc ! Stop time measurement
if (it .gt. nit .and. infinite ) then
    print* , 'ENDGM not converged'
    saveit(method)=it
end if
timeins(method) = time_method(3,method)*60+time_method(4,method)

end subroutine sub_endo