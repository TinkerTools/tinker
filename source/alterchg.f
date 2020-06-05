c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2020 by Chengwen Liu & Jay W. Ponder  ##
c     ##                 All Rights Reserved                  ##
c     ##########################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine alterchg  --  modification of partial charges  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "alterchg" calculates the change in atomic partial charge or
c     monopole values due to bond and angle charge flux coupling
c
c
      subroutine alterchg
      use atoms
      use cflux
      use charge
      use chgpen
      use inform
      use iounit
      use mplpot
      use mpole
      implicit none
      integer i,k
c
c
c     perform dynamic allocation of global array
c
      if (allocated(pcflx)) deallocate (pcflx)
      allocate (pcflx(n))
c
c     zero out the change in charge value at each site
c
      do i = 1, n
         pcflx(i) = 0.0d0
      end do
c
c     find charge modifications due to charge flux
c
      call bndchg
      call angchg
c
c     alter atomic partial charge values for charge flux
c
      if (debug .and. nion.ne.0) then
         write (iout,10)
   10    format (/,' Charge Flux Modification of Partial Charges :',
     &           //,4x,'Atom',14x,'Base Value',6x,'Current',/)
      end if
      do i = 1, nion
         k = iion(i)
         pchg(i) = pchg0(i) + pcflx(k)
         if (debug) then
            write (iout,20)  k,pchg0(i),pchg(i)
   20       format (i8,9x,2f14.5)
         end if
      end do
      
c
c     alter monopoles and charge penetration for charge flux
c
      if (debug .and. npole.ne.0) then
         write (iout,30)
   30    format (/,' Charge Flux Modification of Atomic Monopoles :',
     &           //,4x,'Atom',14x,'Base Value',6x,'Current',/)
      end if
      do i = 1, npole
         k = ipole(i)
         pole(1,i) = mono0(i) + pcflx(k)
         if (use_chgpen)  pval(i) = pval0(i) + pcflx(k)
         if (debug) then
            write (iout,40)  k,mono0(i),pole(1,i)
   40       format (i8,9x,2f14.5)
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine bndchg  --  charge change with bond stretch  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "bndchg" computes modifications to atomic partial charges or
c     monopoles due to bond stretch using a charge flux formulation
c
c
      subroutine bndchg
      use sizes
      use atomid
      use atoms
      use bndstr
      use bound
      use cflux
      use couple
      use mutant
      implicit none
      integer i,j,ia,ib
      integer atoma,atomb
      integer nha,nhb
      integer n12a,n12b
      real*8 xab,yab,zab,rab
      real*8 pjb,pb0,dq
      real*8 priority
      logical muta,mutb
c
c
c     loop over all the bond distances in the system
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         atoma = atomic(ia)
         atomb = atomic(ib)
         pjb = jb(i)
         pb0 = bl(i)
         muta = mut(ia)
         mutb = mut(ib)
         if ((muta) .or. (mutb)) then
            pjb = pjb * elambda
         end if
c
c     determine the higher priority of the bonded atoms
c
         if (atoma .ne. atomb) then
            if (atoma .gt. atomb) then
               priority = 1.0d0
            else
               priority = -1.0d0
            end if
         else
            n12a = n12(ia)
            n12b = n12(ib)
            if (n12a .ne. n12b) then
               if (n12a .gt. n12b) then
                  priority = 1.0d0
               else
                  priority = -1.0d0
               end if
            else
               nha = 0
               nhb = 0
               do j = 1, n12a
                  if (atomic(i12(j,ia)) .eq. 1) then
                     nha = nha + 1
                  end if
               end do
               do j = 1, n12b
                  if (atomic(i12(j,ib)) .eq. 1) then
                     nhb = nhb + 1
                  end if
               end do
               if (nha .ne. nhb) then
                  if (nha .gt. nhb) then
                     priority = 1.0d0
                  else
                     priority = -1.0d0
                  end if
               else
                  priority = 0.0d0
               end if
            end if
         end if
c
c     compute the bond length value for the current bond
c
         xab = x(ia) - x(ib)
         yab = y(ia) - y(ib)
         zab = z(ia) - z(ib)
         if (use_polymer)  call image (xab,yab,zab)
         rab = sqrt(xab*xab + yab*yab + zab*zab)
c
c     find the charge flux increment for the current bond
c
         dq = pjb * (rab-pb0)
         pcflx(ia) = pcflx(ia) - dq*priority
         pcflx(ib) = pcflx(ib) + dq*priority
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine angchg  --  charge change with angle bending  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "angchg" computes modifications to atomic partial charges or
c     monopoles due to angle bending using a charge flux formulation
c
c
      subroutine angchg
      use sizes
      use angbnd
      use atoms
      use bound
      use cflux
      use math
      use mutant
      implicit none
      integer i,ia,ib,ic
      real*8 angle
      real*8 rab,rcb
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 dot,cosine
      real*8 ptheta0
      real*8 pb10,pb20
      real*8 pjbp1,pjbp2
      real*8 pjtheta1
      real*8 pjtheta2
      real*8 dq1,dq2
      logical muta,mutb,mutc
c
c
c     loop over all the bond angles in the system
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
c
c     assign the charge flux parameters for this angle
c
         ptheta0 = anat(i)
         pb10 = bp0(1,i)
         pb20 = bp0(2,i)
         pjbp1 = jbp(1,i)
         pjbp2 = jbp(2,i)
         pjtheta1 = jtheta(1,i)
         pjtheta2 = jtheta(2,i)
         muta = mut(ia)
         mutb = mut(ib)
         mutc = mut(ic)
         if ((muta) .or. (mutb) .or. (mutc)) then
            pjbp1 = pjbp1 * elambda
            pjbp2 = pjbp2 * elambda
            pjtheta1 = pjtheta1 * elambda
            pjtheta2 = pjtheta2 * elambda
         end if
c
c     calculate the angle values and included bond lengths
c
         xia = x(ia)
         yia = y(ia)
         zia = z(ia)
         xib = x(ib)
         yib = y(ib)
         zib = z(ib)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xab = xia - xib
         yab = yia - yib
         zab = zia - zib
         xcb = xic - xib
         ycb = yic - yib
         zcb = zic - zib
         if (use_polymer) then
            call image (xab,yab,zab)
            call image (xcb,ycb,zcb)
         end if
         rab = sqrt(xab*xab + yab*yab + zab*zab)
         rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
         if (rab.ne.0.0d0 .and. rcb.ne.0.0d0) then
            dot = xab*xcb + yab*ycb + zab*zcb
            cosine = dot / (rab*rcb)
            cosine = min(1.0d0,max(-1.0d0,cosine))
            angle = radian * acos(cosine)
         end if
c
c     find the charge flux increment for the current angle
c
         dq1 = pjbp1*(rcb-pb20) + pjtheta1*(angle-ptheta0)
         dq2 = pjbp2*(rab-pb10) + pjtheta2*(angle-ptheta0)
         pcflx(ia) = pcflx(ia) + dq1
         pcflx(ic) = pcflx(ic) + dq2
         pcflx(ib) = pcflx(ib) - dq1 - dq2
      end do
      return
      end
