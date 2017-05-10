!ALBAUGH
      subroutine auxrespa (istep,dt)
      use sizes
      use atomid
      use atoms
      use freeze
      use ielscf
      use moldyn
      use polar
      use units
      use usage
      use virial
      use potent
      implicit none
      integer i,j,k,l
      integer istep
      integer nalt
      integer nmed
      real*8 dt,dt_2
      real*8 dta,dta_2
      real*8 dtm,dtm_2
      real*8 epot,etot
      real*8 eksum,eps
      real*8 temp,pres
      real*8 ealt,dalt
      real*8 dmed
      real*8 term
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 viralt(3,3)
      real*8 ep,es,ef
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: derivsf(:,:)
      real*8, allocatable :: derivss(:,:)
c
c
c     set some time values for the dynamics integration
c
      eps =  0.00000001d0
      dmed = 0.001d0
      nmed = int(dt/(dmed+eps)) + 1
      dmed = dble(nmed)
      dt_2 = 0.5d0 * dt
      dtm = dt / dmed
      dtm_2 = 0.5d0 * dtm
      dalt = 0.00025d0
      nalt = int(dtm/(dalt+eps)) + 1
      dalt = dble(nalt)
      dta = dtm / dalt
      dta_2 = 0.5d0 * dta
      
      allocate(derivs(3,n))
      allocate(derivsf(3,n))
      allocate(derivss(3,n))
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
c
c     store the current atom positions, then find half-step
c     velocities via velocity Verlet recursion
c
      call temper2 (dt,temp)
      
      do l = 1, nmed
         do k = 1, nalt
            do i = 1, n
               if (use(i)) then
                  do j = 1, 3
                     vaux(j,i) = vaux(j,i) + aaux(j,i)*dta_2
                     vpaux(j,i) = vpaux(j,i) + apaux(j,i)*dta_2
                     v(j,i) = v(j,i) + a(j,i)*dta_2
                     uaux(j,i) = uaux(j,i) + vaux(j,i)*dta
                     upaux(j,i) = upaux(j,i) + vpaux(j,i)*dta
                  end do
                  xold(i) = x(i)
                  yold(i) = y(i)
                  zold(i) = z(i)
                  x(i) = x(i) + v(1,i)*dta
                  y(i) = y(i) + v(2,i)*dta
                  z(i) = z(i) + v(3,i)*dta
               end if
            end do
            if (use_rattle)  call rattle (dta,xold,yold,zold)
            call gradfast(ef,derivsf)
            do i = 1, n
               if (use(i)) then
                  do j = 1, 3
                     derivs(j,i) = derivsf(j,i)
                     aaux(j,i) = 0.0d0
                     apaux(j,i) = 0.0d0
                  end do
               end if
            end do
            
            do i = 1, 3
               do j = 1, 3
                  viralt(j,i) = viralt(j,i) + 
     &                          vir(j,i)/(dble(nalt)*dble(nmed))
               end do
            end do
            
            if (k .eq. nalt) then
               call gradint
               term = 2.0d0 / (dtm*dtm)
               if (use_ielscf) then
                  do i = 1, n
                     if (use(i)) then
                        do j = 1, 3
                           aaux(j,i)=term*(uind(j,i)
     &                                   -uaux(j,i))*dble(nalt)
                           apaux(j,i)=term*(uinp(j,i)
     &                                   -upaux(j,i))*dble(nalt)
                        end do
                     end if
                  end do
               else
                  do i = 1, n
                     if (use(i)) then
                        do j = 1, 3
                           aaux(j,i)=gamma_aux*term*(uind(j,i)
     &                                   -uaux(j,i))*dble(nalt)
                           apaux(j,i)=gamma_aux*term*(uinp(j,i)
     &                                   -upaux(j,i))*dble(nalt)
                        end do
                     end if
                  end do
               end if
               if(l .eq. nmed) then
                  call gradslow(es,derivss)
                  do i = 1, n
                     if (use(i)) then
                        do j = 1, 3
                           derivs(j,i)=derivs(j,i)
     &                        +derivss(j,i)*dble(nalt)*dble(nmed)
                        end do
                     end if
                  end do
                  do i = 1, 3
                     do j = 1, 3
                        viralt(j,i) = viralt(j,i) + vir(j,i)
                     end do
                  end do
               end if
            end if
            do i = 1, n
               if (use(i)) then
                  do j = 1, 3
                     a(j,i) = -convert * derivs(j,i) / mass(i)
                     v(j,i) = v(j,i) + a(j,i)*dta_2
                     vaux(j,i) = vaux(j,i) + aaux(j,i)*dta_2
                     vpaux(j,i) = vpaux(j,i) + apaux(j,i)*dta_2
                  end do
               end if
            end do
            if (use_rattle)  call rattle2 (dta)
         end do
      end do
      
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = viralt(j,i)
         end do
      end do
      
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress)
      epot = ef + es
      etot = eksum + epot
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot,eksum)
      call mdrest (istep)
      
      deallocate (derivs)
      deallocate (derivsf)
      deallocate (derivss)
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      return
      end
      
      
      
      subroutine gradint
      implicit none
       
      call chkpole
      call rotpole
      call induce
      return
      end
