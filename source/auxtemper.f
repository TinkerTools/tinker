      subroutine auxtemper (dt)
      use atoms
      use iELSCF
      use units
      implicit none
      integer i,j,nc,ns
      real*8 dt,dtc,dts
      real*8 dt2,dt4,dt8
      real*8 scale,expterm
      real*8 w(3)!w(5)
      
      if(auxstat .eq. 'NOSE-HOOVER') then
c
c     make half-step velocity correction for Nose-Hoover auxiliary dipole system   
c    
         call auxkinetic(.false.)
         nc = 5
         ns = 3
         dtc = dt / dble(nc)
         w(1) = 1.0d0 / (2.0d0-2.0d0**(1.0d0/3.0d0))
         w(2) = 1.0d0 - 2.0d0*w(1)
         w(3) = w(1)
!         w(1) = 1.0d0 / (4.0d0-4.0d0**(1.0d0/3.0d0))
!         w(2) = w(1)
!         w(3) = 1.0d0 - 4.0d0*w(1)
!         w(4) = w(1)
!         w(5) = w(1)
         scale = 1.0d0
         do i = 1, nc
            do j = 1, ns
             dts = w(j) * dtc
             dt2 = 0.5d0 * dts
             dt4 = 0.25d0 * dts
             dt8 = 0.125d0 * dts
            
             gnh_aux(4)=(qnh_aux(3)*vnh_aux(3)*vnh_aux(3)-aux_kelvin)
     &                   /qnh_aux(4)
             vnh_aux(4)=vnh_aux(4)+gnh_aux(4)*dt4
             gnh_aux(3)=(qnh_aux(2)*vnh_aux(2)*vnh_aux(2)-aux_kelvin)
     &                    /qnh_aux(3)
             expterm=exp(-vnh_aux(4)*dt8)
                
             vnh_aux(3)=expterm*(vnh_aux(3)*expterm+gnh_aux(3)*dt4)
             gnh_aux(2)=(qnh_aux(1)*vnh_aux(1)*vnh_aux(1)-aux_kelvin)
     &                   /qnh_aux(2)
             expterm=exp(-vnh_aux(3)*dt8)
               
             vnh_aux(2)=expterm*(vnh_aux(2)*expterm+gnh_aux(2)*dt4)!1/ps
             gnh_aux(1)=(2.0d0*aux_eksum-dble(auxDoF)*aux_kelvin)
     &                   /qnh_aux(1)!1/ps**2
             expterm=exp(-vnh_aux(2)*dt8)
                
             vnh_aux(1)=expterm*(vnh_aux(1)*expterm+gnh_aux(1)*dt4)
             scale=scale*exp(-vnh_aux(1)*dt2)
             aux_eksum=aux_eksum*exp(-vnh_aux(1)*dt2)*
     &                           exp(-vnh_aux(1)*dt2)
                
             pnh_aux(1) = pnh_aux(1) + vnh_aux(1)*dt2!pnhEL dimensionless, vnhEL in 1/ps
             pnh_aux(2) = pnh_aux(2) + vnh_aux(2)*dt2
             pnh_aux(3) = pnh_aux(3) + vnh_aux(3)*dt2
             pnh_aux(4) = pnh_aux(4) + vnh_aux(4)*dt2
               
             gnh_aux(1)=(2.0d0*aux_eksum-dble(auxDoF)*aux_kelvin)
     &                   /qnh_aux(1)
             expterm = exp(-vnh_aux(2)*dt8)
                
             vnh_aux(1)=expterm*(vnh_aux(1)*expterm+gnh_aux(1)*dt4)
             gnh_aux(2)=(qnh_aux(1)*vnh_aux(1)*vnh_aux(1)-aux_kelvin)
     &                   /qnh_aux(2)
             expterm=exp(-vnh_aux(3)*dt8)
                
             vnh_aux(2)=expterm*(vnh_aux(2)*expterm+gnh_aux(2)*dt4)
             gnh_aux(3)=(qnh_aux(2)*vnh_aux(2)*vnh_aux(2)-aux_kelvin)
     &                   /qnh_aux(3)
             expterm=exp(-vnh_aux(4)*dt8)
                
             vnh_aux(3)=expterm*(vnh_aux(3)*expterm+gnh_aux(3)*dt4)
             gnh_aux(4)=(qnh_aux(3)*vnh_aux(3)*vnh_aux(3)-aux_kelvin)
     &                    /qnh_aux(4)
             vnh_aux(4)=vnh_aux(4)+gnh_aux(4)*dt4
            end do
         end do
         
         do i = 1, n
            do j = 1, 3
               v_aux(j,i) = v_aux(j,i)*scale
            end do
         end do
      else if(auxstat .eq. 'NOSE-HOOVER1') then
         call auxkinetic(.false.)
         nc = 5
         ns = 3
         dtc = dt / dble(nc)
         w(1) = 1.0d0 / (2.0d0-2.0d0**(1.0d0/3.0d0))
         w(2) = 1.0d0 - 2.0d0*w(1)
         w(3) = w(1)
!         w(1) = 1.0d0 / (4.0d0-4.0d0**(1.0d0/3.0d0))
!         w(2) = w(1)
!         w(3) = 1.0d0 - 4.0d0*w(1)
!         w(4) = w(1)
!         w(5) = w(1)
         scale = 1.0d0
         do i = 1, nc
            do j = 1, ns
             dts = w(j) * dtc
             dt2 = 0.5d0 * dts
             dt4 = 0.25d0 * dts
             dt8 = 0.125d0 * dts
             
             gnh_aux(1)=(2.0d0*aux_eksum-dble(auxDoF)*aux_kelvin)
     &                    /qnh_aux(1)!1/ps**2
                
             vnh_aux(1)=vnh_aux(1)+gnh_aux(1)*dt4
             scale=scale*exp(-vnh_aux(1)*dt2)
             aux_eksum=aux_eksum*exp(-vnh_aux(1)*dt)
               
             pnh_aux(1) = pnh_aux(1) + vnh_aux(1)*dt2!pnhEL dimensionless, vnhEL in 1/ps
             
             gnh_aux(1)=(2.0d0*aux_eksum-dble(auxDoF)*aux_kelvin)
     &                   /qnh_aux(1)
              
             vnh_aux(1)=(vnh_aux(1)+gnh_aux(1)*dt4)
            end do
         end do
         
         do i = 1, n
            do j = 1, 3
               v_aux(j,i) = v_aux(j,i)*scale
            end do
         end do
      end if
      
      return
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C) 2014 by Alex Albaugh (THG Lab) ##
c     ##                   No Rights Reserved          ##
c     ###################################################
c
c     ############################################################################
c     ##                                                                        ##
c     ##  subroutine dtemper  --  Nose- Hoover thermostat applied at half step  ##
c     ##                                                                        ##
c     ############################################################################
c
c
c     "dtemper" applies a velocity correction at the full time step
c     as needed for the Nose-Hoover extended system thermostat for 
c     real atomic degrees of freedom as well as Drude oscillators
c 
      subroutine auxtemper2 (dt,istep)
      use atoms
      use units
      use iELSCF
      implicit none
      integer i,j,nc,ns
      integer istep
      real*8 scale,dt
      real*8 dtc,dts
      real*8 dt2,dt4,dt8
      real*8 expterm
      real*8 w(3)!w(5)
      
      call auxkinetic(.false.)
      
      if(auxstat .eq. 'NOSE-HOOVER') then
c
c     make full-step velocity correction for Nose-Hoover auxiliary dipole system   
c    
         nc = 5
         ns = 3
         dtc = dt / dble(nc)
         w(1) = 1.0d0 / (2.0d0-2.0d0**(1.0d0/3.0d0))
         w(2) = 1.0d0 - 2.0d0*w(1)
         w(3) = w(1)
!         w(1) = 1.0d0 / (4.0d0-4.0d0**(1.0d0/3.0d0))
!         w(2) = w(1)
!         w(3) = 1.0d0 - 4.0d0*w(1)
!         w(4) = w(1)
!         w(5) = w(1)
         scale = 1.0d0
         do i = 1, nc
            do j = 1, ns
             dts = w(j) * dtc
             dt2 = 0.5d0 * dts
             dt4 = 0.25d0 * dts
             dt8 = 0.125d0 * dts
            
             gnh_aux(4)=(qnh_aux(3)*vnh_aux(3)*vnh_aux(3)-aux_kelvin)
     &                   /qnh_aux(4)
             vnh_aux(4)=vnh_aux(4)+gnh_aux(4)*dt4
             gnh_aux(3)=(qnh_aux(2)*vnh_aux(2)*vnh_aux(2)-aux_kelvin)
     &                   /qnh_aux(3)
             expterm=exp(-vnh_aux(4)*dt8)
                
             vnh_aux(3)=expterm*(vnh_aux(3)*expterm+gnh_aux(3)*dt4)
             gnh_aux(2)=(qnh_aux(1)*vnh_aux(1)*vnh_aux(1)-aux_kelvin)
     &                   /qnh_aux(2)
             expterm=exp(-vnh_aux(3)*dt8)
               
             vnh_aux(2)=expterm*(vnh_aux(2)*expterm+gnh_aux(2)*dt4)!1/ps
             gnh_aux(1)=(2.0d0*aux_eksum-dble(auxDoF)*aux_kelvin)
     &                   /qnh_aux(1)!1/ps**2
             expterm=exp(-vnh_aux(2)*dt8)
                
             vnh_aux(1)=expterm*(vnh_aux(1)*expterm+gnh_aux(1)*dt4)
             scale=scale*exp(-vnh_aux(1)*dt2)
             aux_eksum=aux_eksum*exp(-vnh_aux(1)*dt2)*
     &                           exp(-vnh_aux(1)*dt2)
                
             pnh_aux(1) = pnh_aux(1) + vnh_aux(1)*dt2!pnhEL dimensionless, vnhEL in 1/ps
             pnh_aux(2) = pnh_aux(2) + vnh_aux(2)*dt2
             pnh_aux(3) = pnh_aux(3) + vnh_aux(3)*dt2
             pnh_aux(4) = pnh_aux(4) + vnh_aux(4)*dt2
               
             gnh_aux(1)=(2.0d0*aux_eksum-dble(auxDoF)*aux_kelvin)
     &                   /qnh_aux(1)
             expterm = exp(-vnh_aux(2)*dt8)
                
             vnh_aux(1)=expterm*(vnh_aux(1)*expterm+gnh_aux(1)*dt4)
             gnh_aux(2)=(qnh_aux(1)*vnh_aux(1)*vnh_aux(1)-aux_kelvin)
     &                   /qnh_aux(2)
             expterm=exp(-vnh_aux(3)*dt8)
                
             vnh_aux(2)=expterm*(vnh_aux(2)*expterm+gnh_aux(2)*dt4)
             gnh_aux(3)=(qnh_aux(2)*vnh_aux(2)*vnh_aux(2)-aux_kelvin)
     &                   /qnh_aux(3)
             expterm=exp(-vnh_aux(4)*dt8)
                
             vnh_aux(3)=expterm*(vnh_aux(3)*expterm+gnh_aux(3)*dt4)
             gnh_aux(4)=(qnh_aux(3)*vnh_aux(3)*vnh_aux(3)-aux_kelvin)
     &                    /qnh_aux(4)
             vnh_aux(4)=vnh_aux(4)+gnh_aux(4)*dt4
            end do
         end do
         
         do i = 1, n
            do j = 1, 3
               v_aux(j,i) = v_aux(j,i)*scale
            end do
         end do
      else if(auxstat .eq. 'NOSE-HOOVER1') then
c
c     make full-step velocity correction for Nose-Hoover auxiliary dipole system   
c    
         nc = 5
         ns = 3
         dtc = dt / dble(nc)
         w(1) = 1.0d0 / (2.0d0-2.0d0**(1.0d0/3.0d0))
         w(2) = 1.0d0 - 2.0d0*w(1)
         w(3) = w(1)
!         w(1) = 1.0d0 / (4.0d0-4.0d0**(1.0d0/3.0d0))
!         w(2) = w(1)
!         w(3) = 1.0d0 - 4.0d0*w(1)
!         w(4) = w(1)
!         w(5) = w(1)
         scale = 1.0d0
         do i = 1, nc
            do j = 1, ns
             dts = w(j) * dtc
             dt2 = 0.5d0 * dts
             dt4 = 0.25d0 * dts
             dt8 = 0.125d0 * dts
             
             gnh_aux(1)=(2.0d0*aux_eksum-dble(auxDoF)*aux_kelvin)
     &                    /qnh_aux(1)!1/ps**2
                
             vnh_aux(1)=vnh_aux(1)+gnh_aux(1)*dt4
             scale=scale*exp(-vnh_aux(1)*dt2)
             aux_eksum=aux_eksum*exp(-vnh_aux(1)*dt)
               
             pnh_aux(1) = pnh_aux(1) + vnh_aux(1)*dt2!pnhEL dimensionless, vnhEL in 1/ps
             
             gnh_aux(1)=(2.0d0*aux_eksum-dble(auxDoF)*aux_kelvin)
     &                   /qnh_aux(1)
              
             vnh_aux(1)=(vnh_aux(1)+gnh_aux(1)*dt4)
            end do
         end do
         
         do i = 1, n
            do j = 1, 3
               v_aux(j,i) = v_aux(j,i)*scale
            end do
         end do
      else if (auxstat .eq. 'BERENDSEN') then
         if (aux_temp .eq. 0.0d0)  aux_temp = 0.1d0
         scale=dsqrt(1.0d0+(dt/aux_tautemp)*(aux_kelvin/aux_temp-1.0d0))
         do i = 1, n
            do j = 1, 3
               v_aux(j,i) = scale * v_aux(j,i)
            end do
         end do
      else if (auxstat .eq. 'RESCALE') then
         if(mod(istep,auxscale) .eq. 0) then
            scale = dsqrt(aux_kelvin/aux_temp)
            do i = 1, n
               do j = 1, 3
                  v_aux(j,i) = scale * v_aux(j,i)
               end do
            end do
         end if
      end if

c
c     recompute kinetic energy and instantaneous temperature
c 
      call auxkinetic(.true.)
      return
      end
