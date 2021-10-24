c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine erepel3  --  Pauli repulsion energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erepel3" calculates the Pauli repulsion energy and partitions
c     the energy among the atoms
c
c     literature reference:
c
c     J. A. Rackers and J. W. Ponder, "Classical Pauli Repulsion:
c     An Anisotropic, Atomic Multipole Model", Journal of Chemical
c     Physics, 150, 084104 (2019)
c
c
      subroutine erepel3
      use limits
      use potent
      use reppot
      implicit none
c
c
c     compute the induced dipoles at each polarizable atom
c
      if (reppolar .and. .not.use_polar) then
         use_polar = .true.
         call induce
         use_polar = .false.
      end if
c
c     choose the method for summing over pairwise interactions
c
      if (use_mlist) then
         call erepel3b
      else
         call erepel3a
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine erepel3a  --  Pauli repulsion analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "erepel3a" calculates the Pauli repulsion energy and also
c     partitions the energy among the atoms using a double loop
c
c
      subroutine erepel3a
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mpole
      use mutant
      use polar
      use potent
      use repel
      use reppot
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 e,eterm
      real*8 fgrp,taper
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4,r5
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr,dik,qik
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 diqk,dkqi,qiqk
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 sizi,sizk,sizik
      real*8 vali,valk
      real*8 dmpi,dmpk
      real*8 dmpik(9)
      real*8, allocatable :: rscale(:)
      logical proceed,usei
      logical muti,mutk,mutik
      logical header,huge
      character*6 mode
c
c
c     zero out the repulsion energy and partitioning terms
c
      ner = 0
      er = 0.0d0
      do i = 1, n
         aer(i) = 0.0d0
      end do
      if (nrep .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     perform dynamic allocation of some local arrays
c
      allocate (rscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         rscale(i) = 1.0d0
      end do
c
c     set cutoff distances and switching coefficients
c
      mode = 'REPULS'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug.ne.0 .and. nrep.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Pauli Repulsion Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     calculate the Pauli repulsion interaction energy term
c
      do ii = 1, npole-1
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         sizi = sizpr(ii)
         dmpi = dmppr(ii)
         vali = elepr(ii)
         ci = rpole(1,ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         if (reppolar) then
            dix = dix + uind(1,ii)
            diy = diy + uind(2,ii)
            diz = diz + uind(3,ii)
         end if
         qixx = rpole(5,ii)
         qixy = rpole(6,ii)
         qixz = rpole(7,ii)
         qiyy = rpole(9,ii)
         qiyz = rpole(10,ii)
         qizz = rpole(13,ii)
         usei = use(i)
         muti = mut(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            rscale(i12(j,i)) = r2scale
         end do
         do j = 1, n13(i)
            rscale(i13(j,i)) = r3scale
         end do
         do j = 1, n14(i)
            rscale(i14(j,i)) = r4scale
         end do
         do j = 1, n15(i)
            rscale(i15(j,i)) = r5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, npole
            k = ipole(kk)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  sizk = sizpr(kk)
                  dmpk = dmppr(kk)
                  valk = elepr(kk)
                  ck = rpole(1,kk)
                  dkx = rpole(2,kk)
                  dky = rpole(3,kk)
                  dkz = rpole(4,kk)
                  if (reppolar) then
                     dkx = dkx + uind(1,kk)
                     dky = dky + uind(2,kk)
                     dkz = dkz + uind(3,kk)
                  end if
                  qkxx = rpole(5,kk)
                  qkxy = rpole(6,kk)
                  qkxz = rpole(7,kk)
                  qkyy = rpole(9,kk)
                  qkyz = rpole(10,kk)
                  qkzz = rpole(13,kk)
c
c     intermediates involving moments and separation distance
c
                  dir = dix*xr + diy*yr + diz*zr
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qir = qix*xr + qiy*yr + qiz*zr
                  dkr = dkx*xr + dky*yr + dkz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qkr = qkx*xr + qky*yr + qkz*zr
                  dik = dix*dkx + diy*dky + diz*dkz
                  qik = qix*qkx + qiy*qky + qiz*qkz
                  diqk = dix*qkx + diy*qky + diz*qkz
                  dkqi = dkx*qix + dky*qiy + dkz*qiz
                  qiqk = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                      + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     get reciprocal distance terms for this interaction
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
c
c     get damping coefficients for the Pauli repulsion energy
c
                  call damprep (r,r2,rr1,rr3,rr5,rr7,rr9,rr11,
     &                             9,dmpi,dmpk,dmpik) 
c
c     compute the Pauli repulsion energy for this interaction
c
                  term1 = vali*valk
                  term2 = valk*dir - vali*dkr + dik
                  term3 = vali*qkr + valk*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  eterm = term1*dmpik(1) + term2*dmpik(3)
     &                       + term3*dmpik(5) + term4*dmpik(7)
     &                       + term5*dmpik(9)
                  sizik = sizi * sizk
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = .false.
                  if (muti .or. mutk) then
                     if (vcouple .eq. 1) then
                        mutik = .true.
                     else if (.not.muti .or. .not.mutk) then
                        mutik = .true.
                     end if
                  end if
c
c     get interaction energy, via soft core lambda scaling as needed
c
                  if (mutik) then
                     e = vlambda * sizik * rscale(k) * eterm
     &                      / sqrt(1.0d0-vlambda+r2)
                  else
                     e = sizik * rscale(k) * eterm * rr1
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall Pauli repulsion energy component
c
                  if (e .ne. 0.0d0) then
                     ner = ner + 1
                     er = er + e
                     aer(i) = aer(i) + 0.5d0*e
                     aer(k) = aer(k) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 20.0d0)
                  if ((debug.ne.0 .and. e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Pauli Repulsion',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),r,e
   30                format (' Repuls',4x,2(i7,'-',a3),
     &                          9x,f10.4,2x,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            rscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            rscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            rscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            rscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction energy with other unit cells
c
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            sizi = sizpr(ii)
            dmpi = dmppr(ii)
            vali = elepr(ii)
            ci = rpole(1,ii)
            dix = rpole(2,ii)
            diy = rpole(3,ii)
            diz = rpole(4,ii)
            if (reppolar) then
               dix = dix + uind(1,ii)
               diy = diy + uind(2,ii)
               diz = diz + uind(3,ii)
            end if
            qixx = rpole(5,ii)
            qixy = rpole(6,ii)
            qixz = rpole(7,ii)
            qiyy = rpole(9,ii)
            qiyz = rpole(10,ii)
            qizz = rpole(13,ii)
            usei = use(i)
            muti = mut(i)
c
c     set exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               rscale(i12(j,i)) = r2scale
            end do
            do j = 1, n13(i)
               rscale(i13(j,i)) = r3scale
            end do
            do j = 1, n14(i)
               rscale(i14(j,i)) = r4scale
            end do
            do j = 1, n15(i)
               rscale(i15(j,i)) = r5scale
            end do
c
c     evaluate all sites within the cutoff distance
c
            do kk = ii, npole
               k = ipole(kk)
               mutk = mut(k)
               proceed = .true.
               if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
               if (.not. use_intra)  proceed = .true.
               if (proceed)  proceed = (usei .or. use(k))
               if (proceed) then
                  do jcell = 2, ncell
                     xr = x(k) - xi
                     yr = y(k) - yi
                     zr = z(k) - zi
                     call imager (xr,yr,zr,jcell)
                     r2 = xr*xr + yr* yr + zr*zr
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        sizk = sizpr(kk)
                        dmpk = dmppr(kk)
                        valk = elepr(kk)
                        ck = rpole(1,kk)
                        dkx = rpole(2,kk)
                        dky = rpole(3,kk)
                        dkz = rpole(4,kk)
                        if (reppolar) then
                           dkx = dkx + uind(1,kk)
                           dky = dky + uind(2,kk)
                           dkz = dkz + uind(3,kk)
                        end if
                        qkxx = rpole(5,kk)
                        qkxy = rpole(6,kk)
                        qkxz = rpole(7,kk)
                        qkyy = rpole(9,kk)
                        qkyz = rpole(10,kk)
                        qkzz = rpole(13,kk)
c
c     intermediates involving moments and separation distance
c
                        dir = dix*xr + diy*yr + diz*zr
                        qix = qixx*xr + qixy*yr + qixz*zr
                        qiy = qixy*xr + qiyy*yr + qiyz*zr
                        qiz = qixz*xr + qiyz*yr + qizz*zr
                        qir = qix*xr + qiy*yr + qiz*zr
                        dkr = dkx*xr + dky*yr + dkz*zr
                        qkx = qkxx*xr + qkxy*yr + qkxz*zr
                        qky = qkxy*xr + qkyy*yr + qkyz*zr
                        qkz = qkxz*xr + qkyz*yr + qkzz*zr
                        qkr = qkx*xr + qky*yr + qkz*zr
                        dik = dix*dkx + diy*dky + diz*dkz
                        qik = qix*qkx + qiy*qky + qiz*qkz
                        diqk = dix*qkx + diy*qky + diz*qkz
                        dkqi = dkx*qix + dky*qiy + dkz*qiz
                        qiqk = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                            + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     get reciprocal distance terms for this interaction
c
                        rr1 = 1.0d0 / r
                        rr3 = rr1 / r2
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr9 = 7.0d0 * rr7 / r2
c
c     get damping coefficients for the Pauli repulsion energy
c
                        call damprep (r,r2,rr1,rr3,rr5,rr7,rr9,rr11,
     &                                   9,dmpi,dmpk,dmpik)                  
c
c     compute the Pauli repulsion energy for this interaction
c
                        term1 = vali*valk
                        term2 = valk*dir - vali*dkr + dik
                        term3 = vali*qkr + valk*qir - dir*dkr
     &                             + 2.0d0*(dkqi-diqk+qiqk)
                        term4 = dir*qkr - dkr*qir - 4.0d0*qik
                        term5 = qir*qkr
                        eterm = term1*dmpik(1) + term2*dmpik(3)
     &                             + term3*dmpik(5) + term4*dmpik(7)
     &                             + term5*dmpik(9)
                        sizik = sizi * sizk
c
c     set use of lambda scaling for decoupling or annihilation
c
                        mutik = .false.
                        if (muti .or. mutk) then
                           if (vcouple .eq. 1) then
                              mutik = .true.
                           else if (.not.muti .or. .not.mutk) then
                              mutik = .true.
                           end if
                        end if
c
c     get interaction energy, via soft core lambda scaling as needed
c
                        if (mutik) then
                           e = vlambda * sizik * rscale(k) * eterm
     &                            / sqrt(1.0d0-vlambda+r2)
                        else
                           e = sizik * rscale(k) * eterm * rr1
                        end if
c
c     use energy switching if near the cutoff distance
c
                        if (r2 .gt. cut2) then
                           r3 = r2 * r
                           r4 = r2 * r2
                           r5 = r2 * r3
                           taper = c5*r5 + c4*r4 + c3*r3
     &                                + c2*r2 + c1*r + c0
                           e = e * taper
                        end if
c
c     scale the interaction based on its group membership
c
                        if (use_group)  e = e * fgrp
c
c     increment the overall Pauli repulsion energy component;
c     interaction of an atom with its own image counts half
c
                        if (e .ne. 0.0d0) then
                           ner = ner + 1
                           if (i .eq. k) then
                              er = er + 0.5d0*e
                              aer(i) = aer(i) + 0.5d0*e
                           else
                              er = er + e
                              aer(i) = aer(i) + 0.5d0*e
                              aer(k) = aer(k) + 0.5d0*e
                           end if
                        end if
c
c     increment the total intermolecular energy
c
                        einter = einter + e
c
c     print a message if the energy of this interaction is large
c
                        huge = (e .gt. 20.0d0)
                        if ((debug.ne.0 .and. e.ne.0.0d0)
     &                        .or. (verbose.and.huge)) then
                           if (header) then
                              header = .false.
                              write (iout,40)
   40                         format (/,' Individual Pauli Repulsion',
     &                                   ' Interactions :',
     &                                //,' Type',14x,'Atom Names',
     &                                   15x,'Distance',8x,'Energy',/)
                           end if
                           write (iout,50)  i,name(i),k,name(k),r,e
   50                      format (' Repuls',4x,2(i7,'-',a3),1x,
     &                                '(XTAL)',2x,f10.4,2x,f12.4)
                        end if
                     end if
                  end do
               end if
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               rscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               rscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               rscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               rscale(i15(j,i)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (rscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine erepel3b  --  Pauli repulsion analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "erepel3b" calculates the Pauli repulsion energy and also
c     partitions the energy among the atoms using a neighbor list
c
c
      subroutine erepel3b
      use action
      use analyz
      use atoms
      use atomid
      use bound
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mpole
      use mutant
      use neigh
      use polar
      use repel
      use reppot
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 e,eterm
      real*8 fgrp,taper
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4,r5
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr,dik,qik
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 diqk,dkqi,qiqk
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 sizi,sizk,sizik
      real*8 vali,valk
      real*8 dmpi,dmpk
      real*8 dmpik(9)
      real*8, allocatable :: rscale(:)
      logical proceed,usei
      logical muti,mutk,mutik
      logical header,huge
      character*6 mode
c
c
c     zero out Pauli repulsion energy and partitioning terms
c
      ner = 0
      er = 0.0d0
      do i = 1, n
         aer(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     perform dynamic allocation of some local arrays
c
      allocate (rscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         rscale(i) = 1.0d0
      end do
c
c     set cutoff distances and switching coefficients
c
      mode = 'REPULS'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug.ne.0 .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Pauli Repulsion Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,sizpr,dmppr,elepr,rpole,uind,n12,
!$OMP& i12,n13,i13,n14,i14,n15,i15,r2scale,r3scale,r4scale,r5scale,
!$OMP& nelst,elst,use,use_group,use_intra,use_bounds,reppolar,
!$OMP& vcouple,vlambda,mut,cut2,off2,c0,c1,c2,c3,c4,c5,molcule,
!$OMP& name,debug,verbose,header,iout)
!$OMP& firstprivate(rscale)
!$OMP& shared (er,ner,aer,einter)
!$OMP DO reduction(+:er,ner,aer,einter) schedule(guided)
c
c     calculate the Pauli repulsion interaction energy term
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         sizi = sizpr(ii)
         dmpi = dmppr(ii)
         vali = elepr(ii)
         ci = rpole(1,ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         if (reppolar) then
            dix = dix + uind(1,ii)
            diy = diy + uind(2,ii)
            diz = diz + uind(3,ii)
         end if
         qixx = rpole(5,ii)
         qixy = rpole(6,ii)
         qixz = rpole(7,ii)
         qiyy = rpole(9,ii)
         qiyz = rpole(10,ii)
         qizz = rpole(13,ii)
         usei = use(i)
         muti = mut(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            rscale(i12(j,i)) = r2scale
         end do
         do j = 1, n13(i)
            rscale(i13(j,i)) = r3scale
         end do
         do j = 1, n14(i)
            rscale(i14(j,i)) = r4scale
         end do
         do j = 1, n15(i)
            rscale(i15(j,i)) = r5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = ipole(kk)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  sizk = sizpr(kk)
                  dmpk = dmppr(kk)
                  valk = elepr(kk)
                  ck = rpole(1,kk)
                  dkx = rpole(2,kk)
                  dky = rpole(3,kk)
                  dkz = rpole(4,kk)
                  if (reppolar) then
                     dkx = dkx + uind(1,kk)
                     dky = dky + uind(2,kk)
                     dkz = dkz + uind(3,kk)
                  end if
                  qkxx = rpole(5,kk)
                  qkxy = rpole(6,kk)
                  qkxz = rpole(7,kk)
                  qkyy = rpole(9,kk)
                  qkyz = rpole(10,kk)
                  qkzz = rpole(13,kk)
c
c     intermediates involving moments and separation distance
c
                  dir = dix*xr + diy*yr + diz*zr
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qir = qix*xr + qiy*yr + qiz*zr
                  dkr = dkx*xr + dky*yr + dkz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qkr = qkx*xr + qky*yr + qkz*zr
                  dik = dix*dkx + diy*dky + diz*dkz
                  qik = qix*qkx + qiy*qky + qiz*qkz
                  diqk = dix*qkx + diy*qky + diz*qkz
                  dkqi = dkx*qix + dky*qiy + dkz*qiz
                  qiqk = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                      + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     get reciprocal distance terms for this interaction
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
c
c     get damping coefficients for the Pauli repulsion energy
c     
                  call damprep (r,r2,rr1,rr3,rr5,rr7,rr9,rr11,
     &                             9,dmpi,dmpk,dmpik)                  
c
c     compute the Pauli repulsion energy for this interaction
c
                  term1 = vali*valk
                  term2 = valk*dir - vali*dkr + dik
                  term3 = vali*qkr + valk*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  eterm = term1*dmpik(1) + term2*dmpik(3)
     &                       + term3*dmpik(5) + term4*dmpik(7)
     &                       + term5*dmpik(9)
                  sizik = sizi * sizk
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = .false.
                  if (muti .or. mutk) then
                     if (vcouple .eq. 1) then
                        mutik = .true.
                     else if (.not.muti .or. .not.mutk) then
                        mutik = .true.
                     end if
                  end if
c
c     get interaction energy, via soft core lambda scaling as needed
c
                  if (mutik) then
                     e = vlambda * sizik * rscale(k) * eterm
     &                      / sqrt(1.0d0-vlambda+r2)
                  else
                     e = sizik * rscale(k) * eterm * rr1
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall Pauli repulsion energy component
c
                  if (e .ne. 0.0d0) then
                     ner = ner + 1
                     er = er + e
                     aer(i) = aer(i) + 0.5d0*e
                     aer(k) = aer(k) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 20.0d0)
                  if ((debug.ne.0 .and. e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Pauli Repulsion',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),r,e
   30                format (' Repuls',4x,2(i7,'-',a3),
     &                          9x,f10.4,2x,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            rscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            rscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            rscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            rscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (rscale)
      return
      end
