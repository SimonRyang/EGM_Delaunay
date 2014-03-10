subroutine sub_sim

use comtypes
use comparam
use comvar
use comfunc
use comfuncform
use clock

integer none

real(prec),dimension(nsim)::gridsimx, gridsimh

real(prec):: vx,vh
real(prec),dimension(age)::xex,hex,cex,kex,vex,xen,hen,cen,ken,ven,xhy,hhy,chy,khy,vhy

real(prec),dimension(age-1,nsim,nsim)::ecex,ekex,ecen,eken,echy,ekhy


integer sxc,shc,gc

gridsimx(1:nsim)= func_makegrid(x0min,x0max,nsim,curv2)
gridsimh(1:nsim)= func_makegrid(h0min,h0max,nsim,curv2)

ecex=0.0
ekex=0.0
ecen=0.0
eken=0.0
ecey=0.0
ekey=0.0



!print*, '*************** Exogen Method ************* '
!call tic
do sxc= 1, nsim
    do shc = 1,nsim  
        ! EXGM   
        if (algorithm .eq. 1 .or. algorithm .eq. 4) then   
            xex = 0.0
            hex = 0.0
            cex = 0.0
            kex = 0.0
            vex = 0.0

            ! Values of human capital and financial assets in first period        
            xex(1) = gridsimx(sxc)
            hex(1) = gridsimh(shc)
            
            ! For infinite horizon policy functions are the same for all periods 
            calt(:,:) = csimex(:,:,1)
            kalt(:,:) = ksimex(:,:,1)
            valt(:,:) = vsimex(:,:,1)         
            
            ! Computing consumption, investment, financial assets and human capital
            do it = age-1,1,-1
                if (.not. infinite) then    ! Policy functions depend on age of household in finite horizon
                    calt(:,:) = csimex(:,:,it)
                    kalt(:,:) = ksimex(:,:,it)
                    valt(:,:) = vsimex(:,:,it) 
                end if
                
                call sub_exopolation(xex(age-it),hex(age-it),cex(age-it),kex(age-it),vx,vh,vex(age-it))
            
                xex(age-it+1) = (1+ret)*(xex(age-it)+netw*hex((age-it))-cex((age-it))-kex((age-it)))
                hex(age-it+1) = (1-delta)*(hex(age-it)+func_F(kex((age-it))))        
            end do
            
            ! Set last period investment in finite horizon to zero
            if (infinite) then
                call sub_exopolation(xex(age),hex(age),cex(age),kex(age),vx,vh,vex(age))
            else
                cex(age) = xex(age)+netw*hex(age)
                kex(age) = 0.0
                vex(age) = func_U(cex(age))
            end if
            
            ! Euler Error Evaluation
            do it = 1,age-1
                ecex(it,sxc,shc) = func_eulererr_c(cex(it),cex(it+1),hex(it+1))
                ekex(it,sxc,shc) = func_eulererr_k(cex(it),kex(it),hex(it+1),cex(it+1),kex(it+1),vex(it+1)) 
            end do
        end if
    end do
end do
!call toc

!print*, '*************** Endogenous Method ************* '
!call tic
do sxc= 1, nsim
    do shc = 1,nsim         
        ! ENDGM  
        if (algorithm .eq. 2 .or. algorithm .eq. 4) then
            xen = 0.0
            hen = 0.0
            cen = 0.0
            ken = 0.0
            ven = 0.0
            
            ! Values of human capital and financial assets in first period
            xen(1) = gridsimx(sxc)
            hen(1) = gridsimh(shc)
            
            ! For infinite horizon policy functions are the same for all periods
            if (shc .eq. 1 .and. sxc.eq.1 .and. infinite) then
                call triangulation
            end if 
            
            ! Computing consumption, investment, financial assets and human capital
            do it = age-1,1,-1                 
                if (.not. infinite) then ! Policy functions depend on age of household in finite horizon                 
                    if (it.eq.age-1) then
                        call triangulation 
                    else 
                        grid(1,:) = sav_grid(1,:,it+1)
                        grid(2,:) = sav_grid(2,:,it+1)
                        cf(:) = sav_cf(:,it+1)  
                        kf(:) = sav_kf(:,it+1)  
                        vf(:) = sav_vf(:,it+1)              

                        ntri = sav_ntri(it+1)                               
                        ind(:) = sav_ind(:,it+1)                         
                        til(:,:) = sav_til(:,:,it+1)                 
                        tnbr(:,:) = sav_tnbr(:,:,it+1)                    
                    end if
                end if           

                call delaunay_int (xen(age-it),hen(age-it),cen(age-it),ken(age-it),vx,vh,ven(age-it))

                xen(age-it+1) = (1+ret)*(xen(age-it)+netw*hen((age-it))-cen((age-it))-ken((age-it)))
                hen(age-it+1) = (1-delta)*(hen(age-it)+func_F(ken((age-it))))
            end do
            
            ! Set last period investment in finite horizon to zero
            if (infinite) then
                call delaunay_int (xen(age),hen(age),cen(age),ken(age),vx,vh,ven(age))
            else
                cen(age) = xen(age)+netw*hen(age)
                ken(age) = 0.0
                ven(age) = func_U(cen(age))
            end if
            
            ! Euler Error Evaluation
            do it = 1,age-1
                ecen(it,sxc,shc) = func_eulererr_c(cen(it),cen(it+1),hen(it+1))
                eken(it,sxc,shc) = func_eulererr_k(cen(it),ken(it),hen(it+1),cen(it+1),ken(it+1),ven(it+1)) 
            end do
        end if
    end do
end do
!call toc

!print*, '*************** Hybrid Method ************* '       
!call tic
do sxc= 1, nsim
    do shc = 1,nsim         
        ! HEGM
        if (algorithm .eq. 3 .or. algorithm .eq. 4) then  
            xhy = 0.0
            hhy = 0.0
            chy = 0.0
            khy = 0.0
            vhy = 0.0
            
            ! Values of human capital and financial assets in first period
            xhy(1) = gridsimx(sxc)
            hhy(1) = gridsimh(shc)

            ! For infinite horizon policy functions are the same for all periods 
            calt(:,:) = csimhy(:,:,1)
            kalt(:,:) = ksimhy(:,:,1) 
            gridxhy(:,:) = xsimhy(:,:,1)
            valt(:,:) = vsimhy(:,:,1)

            ! Computing consumption, investment, financial assets and human capital       
            do it = age-1,1,-1
                if (.not. infinite) then ! Policy functions depend on age of household in finite horizon
                    calt(:,:) = csimhy(:,:,it)
                    kalt(:,:) = ksimhy(:,:,it) 
                    gridxhy(:,:) = xsimhy(:,:,it)
                    valt(:,:) = vsimhy(:,:,it)
                end if
                
                call sub_hypolation(xhy(age-it),hhy(age-it),chy(age-it),khy(age-it),vx,vh,vhy(age-it))

                xhy(age-it+1) = (1+ret)*(xhy(age-it)+netw*hhy((age-it))-chy((age-it))-khy((age-it)))
                hhy(age-it+1) = (1-delta)*(hhy(age-it)+func_F(khy((age-it))))
            end do

            ! Set last period investment in finite horizon to zero
            if (infinite) then
                call sub_hypolation(xhy(age),hhy(age),chy(age),khy(age),vx,vh,vhy(age))
            else
                chy(age) = xhy(age)+netw*hhy(age)
                khy(age) = 0.0
                vhy(age) = func_U(chy(age))
            end if
            
            ! Euler Error Evaluation        
            do it = 1,age-1
                echy(it,sxc,shc) = func_eulererr_c(chy(it),chy(it+1),hhy(it+1))
                ekhy(it,sxc,shc) = func_eulererr_k(chy(it),khy(it),hhy(it+1),chy(it+1),khy(it+1),vhy(it+1)) 
            end do
        end if
    end do
end do
!call toc
! Maximum Euler Error

mcx = maxval(abs(ecex(:age-1,:,:)))
mkx = maxval(abs(ekex(:age-1,:,:)))
mcn = maxval(abs(ecen(:age-1,:,:)))
mkn = maxval(abs(eken(:age-1,:,:)))
mcy = maxval(abs(echy(:age-1,:,:)))
mky = maxval(abs(ekhy(:age-1,:,:)))

! Average Euler Error
acx = sum(abs(ecex(:age-1,:,:)))/((age-1)*nsim*nsim)
akx = sum(abs(ekex(:age-1,:,:)))/((age-1)*nsim*nsim)
acn = sum(abs(ecen(:age-1,:,:)))/((age-1)*nsim*nsim)
akn = sum(abs(eken(:age-1,:,:)))/((age-1)*nsim*nsim)
acy = sum(abs(echy(:age-1,:,:)))/((age-1)*nsim*nsim)
aky = sum(abs(ekhy(:age-1,:,:)))/((age-1)*nsim*nsim)

! Writing Euler error in txt file
open (unit = 101, file = 'simerr.txt', status = 'old')!,position='append'
write (101, '(12f30.15)') mcx,mcn,mcy,mkx,mkn,mky,acx,acn,acy,akx,akn,aky 
close (101)

! Writing life-cycle profiles in txt files
open (unit = 102, file = 'simx.txt', status = 'old')!,position='append'
do gc = 1,age
    write (102, '(3f30.15)') xex(gc), xen(gc), xhy(gc)
end do
close (102)

open (unit = 103, file = 'simh.txt', status = 'old')!,position='append'
do gc = 1,age
    write(103, '(3f30.15)') hex(gc), hen(gc), hhy(gc)
end do
close (103)

open (unit = 104, file = 'simk.txt', status = 'old')!,position='append'
do gc = 1,age
    write (104, '(3f30.15)') kex(gc), ken(gc), khy(gc)
end do
close (104)

open (unit = 105, file = 'simc.txt', status = 'old')!,position='append'
do gc = 1,age
    write (105, '(3f30.15)') cex(gc), cen(gc), chy(gc)
end do
close (105)

open (unit = 106, file = 'simv.txt', status = 'old')!,position='append'
do gc = 1,age
    write (106, '(3f30.15)') vex(gc), ven(gc), vhy(gc)
end do
close (106)

!!  Writing human capital evolution for different initial values
!open (unit = 107, file = 'simtest.txt', status = 'old')!,position='append'
!do gc=1, age
!    write (107, '(11f30.15)') htest1(gc,1), htest1(gc,2), htest1(gc,3), htest1(gc,4), htest1(gc,5), htest1(gc,6), htest1(gc,7), htest1(gc,8), htest1(gc,9), htest1(gc,10), htest1(gc,11)
!end do
!close (107)

mcx = log(mcx)/log(10.0)
mcn = log(mcn)/log(10.0)
mcy = log(mcy)/log(10.0)
mkx = log(mkx)/log(10.0)
mkn = log(mkn)/log(10.0)
mky = log(mky)/log(10.0)
acx = log(acx)/log(10.0)
acn = log(acn)/log(10.0)
acy = log(acy)/log(10.0)
akx = log(akx)/log(10.0)
akn = log(akn)/log(10.0)
aky = log(aky)/log(10.0)

!print* , '     EXGM     ENDGM      HEGM'
!print* , 'Maximum log10 Euler Error for consumption'
!write (*, '(3f10.2)') mcx,mcn,mcy
!print* , 'Average log10 Euler Error for consumption'
!write (*, '(3f10.2)') acx,acn,acy 
!print* , 'Maximum log10 Euler Error for investment'
!write (*, '(3f10.2)') mkx,mkn,mky 
!print* , 'Average log10 Euler Error for investment'
!write (*, '(3f10.2)') akx,akn,aky 

print* , '     Max_c  Max_i  Average_c Average_i '
print* , 'ENDGM'
write (*, '(4f10.2)') mcn, mkn, acn, akn
print* , 'HEGM'
write (*, '(4f10.2)') mcy, mky, acy, aky 
print* , 'EXGM'
write (*, '(4f10.2)') mcx, mkx, acx, akx

end subroutine sub_sim