c
c
c     #######################################################
c     ##  COPYRIGHT (C)  2014  by  Alex Albaugh (THG Lab)  ##
c     ##              All Rights Reserved                  ##
c     #######################################################
c
c     ######################################################################################
c     ##                                                                                  ##
c     ##  subroutine dkinetic -- compute kinetic energy components for Drude calculations ##
c     ##                                                                                  ##
c     ######################################################################################
c
c
c     "dkinetic" computes the total kinetic energy and kinetic energy
c     contributions to the pressure tensor by summing over velocities
c     This calculation is specifically for Drude polarization as it
c     further divides kinetic energy and components amongst the various
c     Drude degrees of freedom.
c
c
      subroutine auxkinetic(flag)
      use atoms
      use units
      use iELSCF
      use mpole
      implicit none
      logical flag
      integer i,j,k,c
      real*8 auxterm
      real*8 v1,v2
      real*8 v1p,v2p
c
c
c     zero out the total kinetic energy and its outer product
c  
      aux_eksum = 0.0d0
      auxp_eksum = 0.0d0
      do i = 1, 3
         do j = 1, 3
            aux_ekin(j,i) = 0.0d0
            auxp_ekin(j,i) = 0.0d0
         end do
      end do
c
c     get the total kinetic energy and tensor for atomic, Drude, and H Drude sites
c 
      do i = 1, n
         auxterm = 0.5d0
         do j = 1, 3
            v1 = v_aux(j,i)
            v1p = vp_aux(j,i)
            do k = 1, 3
               v2 = v_aux(k,i)
               v2p = vp_aux(k,i)
               aux_ekin(k,j) = aux_ekin(k,j) + auxterm * v1 * v2
               auxp_ekin(k,j) = auxp_ekin(k,j) + auxterm*v1p*v2p
c               do c = 1, 4
c                  aux_ekin(k,j) = aux_ekin(k,j) + q_iso1(c,k,j)*
cc     &                 4.0d0*auxterm*v_iso1(c,k,j)*v_iso1(c,k,j)/5.0d0
c                  auxp_ekin(k,j) = auxp_ekin(k,j) + qp_iso1(c,k,j)*
c     &                 4.0d0*auxterm*vp_iso1(c,k,j)*vp_iso1(c,k,j)/5.0d0
c               end do
            end do
         end do
      end do

      aux_eksum = aux_ekin(1,1) + aux_ekin(2,2) + aux_ekin(3,3)
      auxp_eksum = auxp_ekin(1,1)+auxp_ekin(2,2)+auxp_ekin(3,3)
      aux_temp = 2.0d0 * aux_eksum / dble(auxDoF)
      auxp_temp = 2.0d0*auxp_eksum/dble(auxDoF)
      if(flag) write(140,*) aux_eksum, aux_temp, auxp_eksum, auxp_temp
      if ((aux_eksum.lt.0.0d0) .or. (auxp_eksum.lt.0.0d0)) then
         print*, "NEGATIVE AUXILIARY 'KINETIC' ENERGY!"
         stop
      end if
      return
      end
