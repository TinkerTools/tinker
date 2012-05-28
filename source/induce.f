c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine induce  --  evaluate induced dipole moments  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "induce" computes the induced dipole moments at polarizable
c     sites due to direct or mutual polarization
c
c     assumes multipole components have already been rotated into
c     the global coordinate frame; computes induced dipoles based
c     on full system, use of active or inactive atoms is ignored
c
c
      subroutine induce
      implicit none
      include 'sizes.i'
      include 'cutoff.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'polar.i'
      include 'solute.i'
      include 'units.i'
      include 'uprior.i'
      integer i,j,k
      real*8 norm
      logical header
c
c
c     choose the method for computation of induced dipoles
c
      if (solvtyp(1:2) .eq. 'PB') then
         call induce0f
      else if (solvtyp(1:2) .eq. 'GK') then
         call induce0e
      else if (use_ewald) then
         if (use_mlist) then
            call induce0d
         else
            call induce0c
         end if
      else
         if (use_mlist) then
            call induce0b
         else
            call induce0a
         end if
      end if
c
c     update the lists of previous induced dipole values
c
      if (use_aspc) then
         nualt = min(nualt+1,maxualt)
         do i = 1, npole
            do j = 1, 3
               do k = nualt, 2, -1
                  udalt(k,j,i) = udalt(k-1,j,i)
                  upalt(k,j,i) = upalt(k-1,j,i)
               end do
               udalt(1,j,i) = uind(j,i)
               upalt(1,j,i) = uinp(j,i)
            end do
         end do
      end if
c
c     print out a list of the final induced dipole moments
c
      if (debug) then
         header = .true.
         do i = 1, npole
            if (polarity(i) .ne. 0.0d0) then
               if (header) then
                  header = .false.
                  if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
                     write (iout,10)
   10                format (/,' Vacuum Induced Dipole Moments',
     &                          ' (Debyes) :')
                  else
                     write (iout,20)
   20                format (/,' Induced Dipole Moments (Debyes) :')
                  end if
                  if (digits .ge. 8) then
                     write (iout,30)
   30                format (/,4x,'Atom',14x,'X',15x,'Y',15x,'Z',
     &                          15x,'Total',/)
                  else if (digits .ge. 6) then
                     write (iout,40)
   40                format (/,4x,'Atom',14x,'X',13x,'Y',13x,'Z',
     &                          12x,'Total',/)
                  else
                     write (iout,50)
   50                format (/,4x,'Atom',14x,'X',11x,'Y',11x,'Z',
     &                          9x,'Total',/)
                  end if
               end if
               k = ipole(i)
               norm = sqrt(uind(1,i)**2+uind(2,i)**2+uind(3,i)**2)
               if (digits .ge. 8) then
                  write (iout,60)  k,(debye*uind(j,i),j=1,3),debye*norm
   60             format (i8,3x,4f16.8)
               else if (digits .ge. 6) then
                  write (iout,70)  k,(debye*uind(j,i),j=1,3),debye*norm
   70             format (i8,4x,4f14.6)
               else
                  write (iout,80)  k,(debye*uind(j,i),j=1,3),debye*norm
   80             format (i8,5x,4f12.4)
               end if
            end if
         end do
         header = .true.
         if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
            do i = 1, npole
               if (polarity(i) .ne. 0.0d0) then
                  if (header) then
                     header = .false.
                     write (iout,90)
   90                format (/,' SCRF Induced Dipole Moments',
     &                          ' (Debyes) :')
                     if (digits .ge. 8) then
                        write (iout,100)
  100                   format (/,4x,'Atom',14x,'X',15x,'Y',15x,'Z',
     &                             15x,'Total',/)
                     else if (digits .ge. 6) then
                        write (iout,110)
  110                   format (/,4x,'Atom',14x,'X',13x,'Y',13x,'Z',
     &                             12x,'Total',/)
                     else
                        write (iout,120)
  120                   format (/,4x,'Atom',14x,'X',11x,'Y',11x,'Z',
     &                             9x,'Total',/)
                     end if
                  end if
                  k = ipole(i)
                  norm = sqrt(uinds(1,i)**2+uinds(2,i)**2+uinds(3,i)**2)
                  if (digits .ge. 8) then
                     write (iout,130)  k,(debye*uinds(j,i),j=1,3),
     &                                 debye*norm
  130                format (i8,3x,4f16.8)
                  else if (digits .ge. 6) then
                     write (iout,140)  k,(debye*uinds(j,i),j=1,3),
     &                                 debye*norm
  140                format (i8,4x,4f14.6)
                  else
                     write (iout,150)  k,(debye*uinds(j,i),j=1,3),
     &                                 debye*norm
  150                format (i8,5x,4f12.4)
                  end if
               end if
            end do
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine induce0a  --  induced dipoles via double loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "induce0a" computes the induced dipole moments at polarizable
c     sites using a pairwise double loop
c
c
      subroutine induce0a
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'units.i'
      include 'uprior.i'
      integer i,j,k,m
      integer ii,kk
      integer iter,maxiter
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,duir,puir
      real*8 dkr,dukr,pukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 udsum,upsum
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: uold(:,:)
      real*8, allocatable :: uoldp(:,:)
      logical proceed,done
      character*6 mode
c
c
c     zero out the induced dipoles at each site
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     compute the direct induced dipole moment at each atom
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = i+1, npole
            dscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                        scale7 = 1.0d0 - expdamp
     &                              *(1.0d0-damp+0.6d0*damp**2)
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
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
                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)*dscale(kk)
                     field(j,k) = field(j,k) + fkd(j)*dscale(kk)
                     fieldp(j,i) = fieldp(j,i) + fid(j)*pscale(kk)
                     fieldp(j,k) = fieldp(j,k) + fkd(j)*pscale(kk)
                  end do
               end if
            end if
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do i = 1, npole
            ii = ipole(i)
            pdi = pdamp(i)
            pti = thole(i)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            qixx = rpole(5,i)
            qixy = rpole(6,i)
            qixz = rpole(7,i)
            qiyy = rpole(9,i)
            qiyz = rpole(10,i)
            qizz = rpole(13,i)
            do j = i, npole
               dscale(ipole(j)) = 1.0d0
               pscale(ipole(j)) = 1.0d0
            end do
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = p2scale
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = p3scale
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = p4scale
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = p5scale
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = d1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = d2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = d3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = d4scale
            end do
            do k = i, npole
               kk = ipole(k)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                           scale7 = 1.0d0 - expdamp
     &                                 *(1.0d0-damp+0.6d0*damp**2)
                        end if
                     end if
                     rr3 = scale3 / (r*r2)
                     rr5 = 3.0d0 * scale5 / (r*r2*r2)
                     rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
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
                     fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                     fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                     fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                     fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                     fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                     fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
                     do j = 1, 3
                        fip(j) = fid(j)
                        fkp(j) = fkd(j)
                     end do
                     if (use_polymer .and. r2 .le. polycut2) then
                        do j = 1, 3
                           fid(j) = fid(j) * dscale(kk)
                           fip(j) = fip(j) * pscale(kk)
                           fkd(j) = fkd(j) * dscale(kk)
                           fkp(j) = fkp(j) * pscale(kk)
                        end do
                     end if
                     do j = 1, 3
                        field(j,i) = field(j,i) + fid(j)
                        fieldp(j,i) = fieldp(j,i) + fip(j)
                        if (ii .ne. kk) then
                           field(j,k) = field(j,k) + fkd(j)
                           fieldp(j,k) = fieldp(j,k) + fkp(j)
                        end if
                     end do
                  end if
               end do
            end do
         end do
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (udir(3,npole))
      allocate (udirp(3,npole))
      allocate (uold(3,npole))
      allocate (uoldp(3,npole))
c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         eps = 100.0d0
c
c     predicted values for always stable predictor-corrector method
c
         if (use_aspc .and. nualt.eq.maxualt) then
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  do k = 1, nualt
                     udsum = udsum + baspc(k)*udalt(k,j,i)
                     upsum = upsum + baspc(k)*upalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
               end do
            end do
         end if
c
c     compute mutual induced dipole moments by an iterative method
c
         do while (.not. done)
            do i = 1, npole
               do j = 1, 3
                  field(j,i) = 0.0d0
                  fieldp(j,i) = 0.0d0
               end do
            end do
            do i = 1, npole-1
               ii = ipole(i)
               pdi = pdamp(i)
               pti = thole(i)
               duix = uind(1,i)
               duiy = uind(2,i)
               duiz = uind(3,i)
               puix = uinp(1,i)
               puiy = uinp(2,i)
               puiz = uinp(3,i)
               do j = i+1, npole
                  dscale(ipole(j)) = 1.0d0
               end do
               do j = 1, np11(ii)
                  dscale(ip11(j,ii)) = u1scale
               end do
               do j = 1, np12(ii)
                  dscale(ip12(j,ii)) = u2scale
               end do
               do j = 1, np13(ii)
                  dscale(ip13(j,ii)) = u3scale
               end do
               do j = 1, np14(ii)
                  dscale(ip14(j,ii)) = u4scale
               end do
               do k = i+1, npole
                  kk = ipole(k)
                  proceed = .true.
                  if (use_intra)
     &               call groups (proceed,fgrp,ii,kk,0,0,0,0)
                  if (proceed) then
                     xr = x(kk) - x(ii)
                     yr = y(kk) - y(ii)
                     zr = z(kk) - z(ii)
                     call image (xr,yr,zr)
                     r2 = xr*xr + yr* yr + zr*zr
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        dukx = uind(1,k)
                        duky = uind(2,k)
                        dukz = uind(3,k)
                        pukx = uinp(1,k)
                        puky = uinp(2,k)
                        pukz = uinp(3,k)
                        scale3 = dscale(kk)
                        scale5 = dscale(kk)
                        damp = pdi * pdamp(k)
                        if (damp .ne. 0.0d0) then
                           pgamma = min(pti,thole(k))
                           damp = -pgamma * (r/damp)**3
                           if (damp .gt. -50.0d0) then
                              expdamp = exp(damp)
                              scale3 = scale3 * (1.0d0-expdamp)
                              scale5 = scale5 * (1.0d0-expdamp
     &                                              *(1.0d0-damp))
                           end if
                        end if
                        rr3 = scale3 / (r*r2)
                        rr5 = 3.0d0 * scale5 / (r*r2*r2)
                        duir = xr*duix + yr*duiy + zr*duiz
                        dukr = xr*dukx + yr*duky + zr*dukz
                        puir = xr*puix + yr*puiy + zr*puiz
                        pukr = xr*pukx + yr*puky + zr*pukz
                        fid(1) = -rr3*dukx + rr5*dukr*xr
                        fid(2) = -rr3*duky + rr5*dukr*yr
                        fid(3) = -rr3*dukz + rr5*dukr*zr
                        fkd(1) = -rr3*duix + rr5*duir*xr
                        fkd(2) = -rr3*duiy + rr5*duir*yr
                        fkd(3) = -rr3*duiz + rr5*duir*zr
                        fip(1) = -rr3*pukx + rr5*pukr*xr
                        fip(2) = -rr3*puky + rr5*pukr*yr
                        fip(3) = -rr3*pukz + rr5*pukr*zr
                        fkp(1) = -rr3*puix + rr5*puir*xr
                        fkp(2) = -rr3*puiy + rr5*puir*yr
                        fkp(3) = -rr3*puiz + rr5*puir*zr
                        do j = 1, 3
                           field(j,i) = field(j,i) + fid(j)
                           field(j,k) = field(j,k) + fkd(j)
                           fieldp(j,i) = fieldp(j,i) + fip(j)
                           fieldp(j,k) = fieldp(j,k) + fkp(j)
                        end do
                     end if
                  end if
               end do
            end do
c
c     periodic boundary for large cutoffs via replicates method
c
            if (use_replica) then
               do i = 1, npole
                  ii = ipole(i)
                  pdi = pdamp(i)
                  pti = thole(i)
                  duix = uind(1,i)
                  duiy = uind(2,i)
                  duiz = uind(3,i)
                  puix = uinp(1,i)
                  puiy = uinp(2,i)
                  puiz = uinp(3,i)
                  do j = i, npole
                     dscale(ipole(j)) = 1.0d0
                  end do
                  do j = 1, np11(ii)
                     dscale(ip11(j,ii)) = u1scale
                  end do
                  do j = 1, np12(ii)
                     dscale(ip12(j,ii)) = u2scale
                  end do
                  do j = 1, np13(ii)
                     dscale(ip13(j,ii)) = u3scale
                  end do
                  do j = 1, np14(ii)
                     dscale(ip14(j,ii)) = u4scale
                  end do
                  do k = i, npole
                     kk = ipole(k)
                     dukx = uind(1,k)
                     duky = uind(2,k)
                     dukz = uind(3,k)
                     pukx = uinp(1,k)
                     puky = uinp(2,k)
                     pukz = uinp(3,k)
                     proceed = .true.
                     do m = 1, ncell
                        xr = x(kk) - x(ii)
                        yr = y(kk) - y(ii)
                        zr = z(kk) - z(ii)
                        call imager (xr,yr,zr,m)
                        r2 = xr*xr + yr* yr + zr*zr
                        if (r2 .le. off2) then
                           r = sqrt(r2)
                           scale3 = 1.0d0
                           scale5 = 1.0d0
                           damp = pdi * pdamp(k)
                           if (damp .ne. 0.0d0) then
                              pgamma = min(pti,thole(k))
                              damp = -pgamma * (r/damp)**3
                              if (damp .gt. -50.0d0) then
                                 expdamp = exp(damp)
                                 scale3 = 1.0d0 - expdamp
                                 scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                              end if
                           end if
                           rr3 = scale3 / (r*r2)
                           rr5 = 3.0d0 * scale5 / (r*r2*r2)
                           duir = xr*duix + yr*duiy + zr*duiz
                           dukr = xr*dukx + yr*duky + zr*dukz
                           puir = xr*puix + yr*puiy + zr*puiz
                           pukr = xr*pukx + yr*puky + zr*pukz
                           fid(1) = -rr3*dukx + rr5*dukr*xr
                           fid(2) = -rr3*duky + rr5*dukr*yr
                           fid(3) = -rr3*dukz + rr5*dukr*zr
                           fkd(1) = -rr3*duix + rr5*duir*xr
                           fkd(2) = -rr3*duiy + rr5*duir*yr
                           fkd(3) = -rr3*duiz + rr5*duir*zr
                           fip(1) = -rr3*pukx + rr5*pukr*xr
                           fip(2) = -rr3*puky + rr5*pukr*yr
                           fip(3) = -rr3*pukz + rr5*pukr*zr
                           fkp(1) = -rr3*puix + rr5*puir*xr
                           fkp(2) = -rr3*puiy + rr5*puir*yr
                           fkp(3) = -rr3*puiz + rr5*puir*zr
                           if (use_polymer) then
                              if (r2 .le. polycut2) then
                                 do j = 1, 3
                                    fid(j) = fid(j) * dscale(kk)
                                    fkd(j) = fkd(j) * dscale(kk)
                                    fip(j) = fip(j) * dscale(kk)
                                    fkp(j) = fkp(j) * dscale(kk)
                                 end do
                              end if
                           end if
                           do j = 1, 3
                              field(j,i) = field(j,i) + fid(j)
                              fieldp(j,i) = fieldp(j,i) + fip(j)
                              if (ii .ne. kk) then
                                 field(j,k) = field(j,k) + fkd(j)
                                 fieldp(j,k) = fieldp(j,k) + fkp(j)
                              end if
                           end do
                        end if
                     end do
                  end do
               end do
            end if
c
c     check to see if the mutual induced dipoles have converged
c
            iter = iter + 1
            epsold = eps
            epsd = 0.0d0
            epsp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  uold(j,i) = uind(j,i)
                  uoldp(j,i) = uinp(j,i)
                  uind(j,i) = udir(j,i) + polarity(i)*field(j,i)
                  uinp(j,i) = udirp(j,i) + polarity(i)*fieldp(j,i)
                  uind(j,i) = uold(j,i) + polsor*(uind(j,i)-uold(j,i))
                  uinp(j,i) = uoldp(j,i) + polsor*(uinp(j,i)-uoldp(j,i))
                  epsd = epsd + (uind(j,i)-uold(j,i))**2
                  epsp = epsp + (uinp(j,i)-uoldp(j,i))**2
               end do
            end do
            eps = max(epsd,epsp)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. maxiter)  done = .true.
         end do
         if (debug) then
            write (iout,30)  iter,eps
   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
     &                 6x,'RMS Change',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (eps .gt. poleps) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      deallocate (field)
      deallocate (fieldp)
      deallocate (udir)
      deallocate (udirp)
      deallocate (uold)
      deallocate (uoldp)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine induce0b  --  induced dipoles via neighbor list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "induce0b" computes the induced dipole moments at polarizable
c     sites using a neighbor list
c
c
      subroutine induce0b
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'couple.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'units.i'
      include 'uprior.i'
      integer i,j,k
      integer ii,kk,kkk
      integer iter,maxiter
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,duir,puir
      real*8 dkr,dukr,pukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 udsum,upsum
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: uold(:,:)
      real*8, allocatable :: uoldp(:,:)
      logical proceed,done
      character*6 mode
c
c
c     zero out the induced dipoles at each site
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     compute the direct induced dipole moment at each atom
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, nelst(i)
            dscale(ipole(elst(j,i))) = 1.0d0
            pscale(ipole(elst(j,i))) = 1.0d0
         end do
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                        scale7 = 1.0d0 - expdamp
     &                              *(1.0d0-damp+0.6d0*damp**2)
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
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
                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)*dscale(kk)
                     field(j,k) = field(j,k) + fkd(j)*dscale(kk)
                     fieldp(j,i) = fieldp(j,i) + fid(j)*pscale(kk)
                     fieldp(j,k) = fieldp(j,k) + fkd(j)*pscale(kk)
                  end do
               end if
            end if
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (udir(3,npole))
      allocate (udirp(3,npole))
      allocate (uold(3,npole))
      allocate (uoldp(3,npole))
c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         eps = 100.0d0
c
c     predicted values for always stable predictor-corrector method
c
         if (use_aspc .and. nualt.eq.maxualt) then
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  do k = 1, nualt
                     udsum = udsum + baspc(k)*udalt(k,j,i)
                     upsum = upsum + baspc(k)*upalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
               end do
            end do
         end if
c
c     compute mutual induced dipole moments by an iterative method
c
         do while (.not. done)
            do i = 1, npole
               do j = 1, 3
                  field(j,i) = 0.0d0
                  fieldp(j,i) = 0.0d0
               end do
            end do
            do i = 1, npole
               ii = ipole(i)
               pdi = pdamp(i)
               pti = thole(i)
               duix = uind(1,i)
               duiy = uind(2,i)
               duiz = uind(3,i)
               puix = uinp(1,i)
               puiy = uinp(2,i)
               puiz = uinp(3,i)
               do j = 1, nelst(i)
                  dscale(ipole(elst(j,i))) = 1.0d0
               end do
               do j = 1, np11(ii)
                  dscale(ip11(j,ii)) = u1scale
               end do
               do j = 1, np12(ii)
                  dscale(ip12(j,ii)) = u2scale
               end do
               do j = 1, np13(ii)
                  dscale(ip13(j,ii)) = u3scale
               end do
               do j = 1, np14(ii)
                  dscale(ip14(j,ii)) = u4scale
               end do
               do kkk = 1, nelst(i)
                  k = elst(kkk,i)
                  kk = ipole(k)
                  proceed = .true.
                  if (use_intra)
     &               call groups (proceed,fgrp,ii,kk,0,0,0,0)
                  if (proceed) then
                     xr = x(kk) - x(ii)
                     yr = y(kk) - y(ii)
                     zr = z(kk) - z(ii)
                     call image (xr,yr,zr)
                     r2 = xr*xr + yr* yr + zr*zr
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        dukx = uind(1,k)
                        duky = uind(2,k)
                        dukz = uind(3,k)
                        pukx = uinp(1,k)
                        puky = uinp(2,k)
                        pukz = uinp(3,k)
                        scale3 = dscale(kk)
                        scale5 = dscale(kk)
                        damp = pdi * pdamp(k)
                        if (damp .ne. 0.0d0) then
                           pgamma = min(pti,thole(k))
                           damp = -pgamma * (r/damp)**3
                           if (damp .gt. -50.0d0) then
                              expdamp = exp(damp)
                              scale3 = scale3 * (1.0d0-expdamp)
                              scale5 = scale5 * (1.0d0-(1.0d0-damp)
     &                                              *expdamp)
                           end if
                        end if
                        rr3 = scale3 / (r*r2)
                        rr5 = 3.0d0 * scale5 / (r*r2*r2)
                        duir = xr*duix + yr*duiy + zr*duiz
                        dukr = xr*dukx + yr*duky + zr*dukz
                        puir = xr*puix + yr*puiy + zr*puiz
                        pukr = xr*pukx + yr*puky + zr*pukz
                        fid(1) = -rr3*dukx + rr5*dukr*xr
                        fid(2) = -rr3*duky + rr5*dukr*yr
                        fid(3) = -rr3*dukz + rr5*dukr*zr
                        fkd(1) = -rr3*duix + rr5*duir*xr
                        fkd(2) = -rr3*duiy + rr5*duir*yr
                        fkd(3) = -rr3*duiz + rr5*duir*zr
                        fip(1) = -rr3*pukx + rr5*pukr*xr
                        fip(2) = -rr3*puky + rr5*pukr*yr
                        fip(3) = -rr3*pukz + rr5*pukr*zr
                        fkp(1) = -rr3*puix + rr5*puir*xr
                        fkp(2) = -rr3*puiy + rr5*puir*yr
                        fkp(3) = -rr3*puiz + rr5*puir*zr
                        do j = 1, 3
                           field(j,i) = field(j,i) + fid(j)
                           field(j,k) = field(j,k) + fkd(j)
                           fieldp(j,i) = fieldp(j,i) + fip(j)
                           fieldp(j,k) = fieldp(j,k) + fkp(j)
                        end do
                     end if
                  end if
               end do
            end do
c
c     check to see if the mutual induced dipoles have converged
c
            iter = iter + 1
            epsold = eps
            epsd = 0.0d0
            epsp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  uold(j,i) = uind(j,i)
                  uoldp(j,i) = uinp(j,i)
                  uind(j,i) = udir(j,i) + polarity(i)*field(j,i)
                  uinp(j,i) = udirp(j,i) + polarity(i)*fieldp(j,i)
                  uind(j,i) = uold(j,i) + polsor*(uind(j,i)-uold(j,i))
                  uinp(j,i) = uoldp(j,i) + polsor*(uinp(j,i)-uoldp(j,i))
                  epsd = epsd + (uind(j,i)-uold(j,i))**2
                  epsp = epsp + (uinp(j,i)-uoldp(j,i))**2
               end do
            end do
            eps = max(epsd,epsp)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. maxiter)  done = .true.
         end do
         if (debug) then
            write (iout,30)  iter,eps
   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
     &                 6x,'RMS Change',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (eps .gt. poleps) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      deallocate (field)
      deallocate (fieldp)
      deallocate (udir)
      deallocate (udirp)
      deallocate (uold)
      deallocate (uoldp)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine induce0c  --  Ewald induced dipoles via loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "induce0c" computes the induced dipole moments at polarizable
c     sites using particle mesh Ewald summation a double loop
c
c
      subroutine induce0c
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      include 'units.i'
      include 'uprior.i'
      integer i,j,k,ii
      integer iter,maxiter
      real*8 eps,term
      real*8 epsd,epsp
      real*8 epsold
      real*8 udsum,upsum
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: uold(:,:)
      real*8, allocatable :: uoldp(:,:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      logical done
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (udir(3,npole))
      allocate (udirp(3,npole))
      allocate (uold(3,npole))
      allocate (uoldp(3,npole))
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
c
c     zero out the induced dipole and the field at each site
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     get the reciprical space part of the electrostatic field
c
      call udirect1 (field)
c
c     get the real space portion of the electrostatic field
c
      do i = 1, npole
         do j = 1, 3
            fieldp(j,i) = field(j,i)
         end do
      end do
      call udirect2a (field,fieldp)
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npole
         do j = 1, 3
            field(j,i) = field(j,i) + term*rpole(j+1,i)
            fieldp(j,i) = fieldp(j,i) + term*rpole(j+1,i)
         end do
      end do
c
c     compute the cell dipole boundary correction to field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
         end do
         do i = 1, npole
            ii = ipole(i)
            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucell(j)
            end do
         end do
      end if
c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         eps = 100.0d0
c
c     predicted values for always stable predictor-corrector method
c
         if (use_aspc .and. nualt.eq.maxualt) then
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  do k = 1, nualt
                     udsum = udsum + baspc(k)*udalt(k,j,i)
                     upsum = upsum + baspc(k)*upalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
               end do
            end do
         end if
c
c     compute mutual induced dipole moments by an iterative method
c
         do while (.not. done)
            do i = 1, npole
               do j = 1, 3
                  field(j,i) = 0.0d0
                  fieldp(j,i) = 0.0d0
               end do
            end do
            call umutual1 (field,fieldp)
            call umutual2a (field,fieldp)
            term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
            do i = 1, npole
               do j = 1, 3
                  field(j,i) = field(j,i) + term*uind(j,i)
                  fieldp(j,i) = fieldp(j,i) + term*uinp(j,i)
               end do
            end do
            if (boundary .eq. 'VACUUM') then
               do i = 1, 3
                  ucell(i) = 0.0d0
                  ucellp(i) = 0.0d0
               end do
               do i = 1, npole
                  do j = 1, 3
                     ucell(j) = ucell(j) + uind(j,i)
                     ucellp(j) = ucellp(j) + uinp(j,i)
                  end do
               end do
               term = (4.0d0/3.0d0) * pi/volbox
               do i = 1, npole
                  do j = 1, 3
                     field(j,i) = field(j,i) - term*ucell(j)
                     fieldp(j,i) = fieldp(j,i) - term*ucellp(j)
                  end do
               end do
            end if
c
c     check to see if the mutual induced dipoles have converged
c
            iter = iter + 1
            epsold = eps
            epsd = 0.0d0
            epsp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  uold(j,i) = uind(j,i)
                  uoldp(j,i) = uinp(j,i)
                  uind(j,i) = udir(j,i) + polarity(i)*field(j,i)
                  uinp(j,i) = udirp(j,i) + polarity(i)*fieldp(j,i)
                  uind(j,i) = uold(j,i) + polsor*(uind(j,i)-uold(j,i))
                  uinp(j,i) = uoldp(j,i) + polsor*(uinp(j,i)-uoldp(j,i))
                  epsd = epsd + (uind(j,i)-uold(j,i))**2
                  epsp = epsp + (uinp(j,i)-uoldp(j,i))**2
               end do
            end do
            eps = max(epsd,epsp)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. maxiter)  done = .true.
         end do
         if (debug) then
            write (iout,30)  iter,eps
   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
     &                 6x,'RMS Change',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (eps .gt. poleps) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (udir)
      deallocate (udirp)
      deallocate (uold)
      deallocate (uoldp)
      deallocate (field)
      deallocate (fieldp)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine induce0d  --  Ewald induced dipoles via list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "induce0d" computes the induced dipole moments at polarizable
c     sites using particle mesh Ewald summation and a neighbor list
c
c
      subroutine induce0d
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      include 'units.i'
      include 'uprior.i'
      integer i,j,k,ii
      integer iter,maxiter
      real*8 eps,term
      real*8 epsd,epsp
      real*8 epsold
      real*8 udsum,upsum
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: uold(:,:)
      real*8, allocatable :: uoldp(:,:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      logical done
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (udir(3,npole))
      allocate (udirp(3,npole))
      allocate (uold(3,npole))
      allocate (uoldp(3,npole))
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
c
c     zero out the induced dipole and the field at each site
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     get the reciprical space part of the electrostatic field
c
      call udirect1 (field)
c
c     get the real space portion of the electrostatic field
c
      do i = 1, npole
         do j = 1, 3
            fieldp(j,i) = field(j,i)
         end do
      end do
      call udirect2b (field,fieldp)
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npole
         do j = 1, 3
            field(j,i) = field(j,i) + term*rpole(j+1,i)
            fieldp(j,i) = fieldp(j,i) + term*rpole(j+1,i)
         end do
      end do
c
c     compute the cell dipole boundary correction to field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
         end do
         do i = 1, npole
            ii = ipole(i)
            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucell(j)
            end do
         end do
      end if
c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         eps = 100.0d0
c
c     predicted values for always stable predictor-corrector method
c
         if (use_aspc .and. nualt.eq.maxualt) then
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  do k = 1, nualt
                     udsum = udsum + baspc(k)*udalt(k,j,i)
                     upsum = upsum + baspc(k)*upalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
               end do
            end do
         end if
c
c     compute mutual induced dipole moments by an iterative method
c
         do while (.not. done)
            do i = 1, npole
               do j = 1, 3
                  field(j,i) = 0.0d0
                  fieldp(j,i) = 0.0d0
               end do
            end do
            call umutual1 (field,fieldp)
            call umutual2b (field,fieldp)
            term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
            do i = 1, npole
               do j = 1, 3
                  field(j,i) = field(j,i) + term*uind(j,i)
                  fieldp(j,i) = fieldp(j,i) + term*uinp(j,i)
               end do
            end do
            if (boundary .eq. 'VACUUM') then
               do i = 1, 3
                  ucell(i) = 0.0d0
                  ucellp(i) = 0.0d0
               end do
               do i = 1, npole
                  do j = 1, 3
                     ucell(j) = ucell(j) + uind(j,i)
                     ucellp(j) = ucellp(j) + uinp(j,i)
                  end do
               end do
               term = (4.0d0/3.0d0) * pi/volbox
               do i = 1, npole
                  do j = 1, 3
                     field(j,i) = field(j,i) - term*ucell(j)
                     fieldp(j,i) = fieldp(j,i) - term*ucellp(j)
                  end do
               end do
            end if
c
c     check to see if the mutual induced dipoles have converged
c
            iter = iter + 1
            epsold = eps
            epsd = 0.0d0
            epsp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  uold(j,i) = uind(j,i)
                  uoldp(j,i) = uinp(j,i)
                  uind(j,i) = udir(j,i) + polarity(i)*field(j,i)
                  uinp(j,i) = udirp(j,i) + polarity(i)*fieldp(j,i)
                  uind(j,i) = uold(j,i) + polsor*(uind(j,i)-uold(j,i))
                  uinp(j,i) = uoldp(j,i) + polsor*(uinp(j,i)-uoldp(j,i))
                  epsd = epsd + (uind(j,i)-uold(j,i))**2
                  epsp = epsp + (uinp(j,i)-uoldp(j,i))**2
               end do
            end do
            eps = max(epsd,epsp)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. maxiter)  done = .true.
         end do
         if (debug) then
            write (iout,30)  iter,eps
   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
     &                 6x,'RMS Change',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (eps.gt.poleps .and. .not.use_aspc) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (udir)
      deallocate (udirp)
      deallocate (uold)
      deallocate (uoldp)
      deallocate (field)
      deallocate (fieldp)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine udirect1  --  Ewald recip direct induced field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "udirect1" computes the reciprocal space contribution of the
c     permanent atomic multipole moments to the field
c
c
      subroutine udirect1 (field)
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'boxes.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
      include 'pme.i'
      integer i,j,k,ntot
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff,nf1,nf2,nf3
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 field(3,*)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: cphi(:,:)
      real*8, allocatable :: fphi(:,:)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cmp(10,npole))
      allocate (fmp(10,npole))
      allocate (cphi(10,npole))
      allocate (fphi(20,npole))
c
c     copy multipole moments and coordinates to local storage
c
      do i = 1, npole
         cmp(1,i) = rpole(1,i)
         cmp(2,i) = rpole(2,i)
         cmp(3,i) = rpole(3,i)
         cmp(4,i) = rpole(4,i)
         cmp(5,i) = rpole(5,i)
         cmp(6,i) = rpole(9,i)
         cmp(7,i) = rpole(13,i)
         cmp(8,i) = 2.0d0 * rpole(6,i)
         cmp(9,i) = 2.0d0 * rpole(7,i)
         cmp(10,i) = 2.0d0 * rpole(10,i)
      end do
c
c     compute B-spline coefficients and spatial decomposition
c
      call bspline_fill
      call table_fill
c
c     convert Cartesian multipoles to fractional coordinates
c
      call cmp_to_fmp (cmp,fmp)
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_mpole (fmp)
      call fftfront
c
c     make the scalar summation over reciprocal lattice
c
      qfac(1,1,1) = 0.0d0
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      ntot = nfft1 * nfft2 * nfft3
      do i = 1, ntot-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     account for the zeroth grid point for a finite system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         qfac(1,1,1) = expterm
      end if
c
c     complete the transformation of the PME grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get field
c
      call fftback
      call fphi_mpole (fphi)
c
c     convert the field from fractional to Cartesian
c
      call fphi_to_cphi (fphi,cphi)
c
c     increment the field at each multipole site
c
      do i = 1, npole
         field(1,i) = field(1,i) - cphi(2,i)
         field(2,i) = field(2,i) - cphi(3,i)
         field(3,i) = field(3,i) - cphi(4,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cmp)
      deallocate (fmp)
      deallocate (cphi)
      deallocate (fphi)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine udirect2a  --  Ewald real direct field via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "udirect2a" computes the real space contribution of the permanent
c     atomic multipole moments to the field via a double loop
c
c
      subroutine udirect2a (field,fieldp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr,r,r2
      real*8 erfc,bfac,exp2a
      real*8 drr3,drr5,drr7
      real*8 prr3,prr5,prr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qiyy,qizz
      real*8 qixy,qixz,qiyz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkyy,qkzz
      real*8 qkxy,qkxz,qkyz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 ralpha
      real*8 alsq2,alsq2n
      real*8 pdi,pti,pgamma
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 bn(0:3)
      real*8 fim(3),fkm(3)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      character*6 mode
      external erfc
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 3
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     compute the error function scaled and unscaled terms
c
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                     scale7 = 1.0d0 - expdamp
     &                           *(1.0d0-damp+0.6d0*damp**2)
                  end if
               end if
               dsc3 = scale3 * dscale(kk)
               dsc5 = scale5 * dscale(kk)
               dsc7 = scale7 * dscale(kk)
               psc3 = scale3 * pscale(kk)
               psc5 = scale5 * pscale(kk)
               psc7 = scale7 * pscale(kk)
               drr3 = (1.0d0-dsc3) / (r*r2)
               drr5 = 3.0d0 * (1.0d0-dsc5) / (r*r2*r2)
               drr7 = 15.0d0 * (1.0d0-dsc7) / (r*r2*r2*r2)
               prr3 = (1.0d0-psc3) / (r*r2)
               prr5 = 3.0d0 * (1.0d0-psc5) / (r*r2*r2)
               prr7 = 15.0d0 * (1.0d0-psc7) / (r*r2*r2*r2)
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
               fim(1) = -xr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dkx + 2.0d0*bn(2)*qkx
               fim(2) = -yr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dky + 2.0d0*bn(2)*qky
               fim(3) = -zr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dkz + 2.0d0*bn(2)*qkz
               fkm(1) = xr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*dix - 2.0d0*bn(2)*qix
               fkm(2) = yr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*diy - 2.0d0*bn(2)*qiy
               fkm(3) = zr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*diz - 2.0d0*bn(2)*qiz
               fid(1) = -xr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dkx + 2.0d0*drr5*qkx
               fid(2) = -yr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dky + 2.0d0*drr5*qky
               fid(3) = -zr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dkz + 2.0d0*drr5*qkz
               fkd(1) = xr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*dix - 2.0d0*drr5*qix
               fkd(2) = yr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*diy - 2.0d0*drr5*qiy
               fkd(3) = zr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*diz - 2.0d0*drr5*qiz
               fip(1) = -xr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dkx + 2.0d0*prr5*qkx
               fip(2) = -yr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dky + 2.0d0*prr5*qky
               fip(3) = -zr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dkz + 2.0d0*prr5*qkz
               fkp(1) = xr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*dix - 2.0d0*prr5*qix
               fkp(2) = yr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*diy - 2.0d0*prr5*qiy
               fkp(3) = zr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*diz - 2.0d0*prr5*qiz
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  field(j,i) = field(j,i) + fim(j) - fid(j)
                  field(j,k) = field(j,k) + fkm(j) - fkd(j)
                  fieldp(j,i) = fieldp(j,i) + fim(j) - fip(j)
                  fieldp(j,k) = fieldp(j,k) + fkm(j) - fkp(j)
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do i = 1, npole
            ii = ipole(i)
            pdi = pdamp(i)
            pti = thole(i)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            qixx = rpole(5,i)
            qixy = rpole(6,i)
            qixz = rpole(7,i)
            qiyy = rpole(9,i)
            qiyz = rpole(10,i)
            qizz = rpole(13,i)
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = p2scale
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = p3scale
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = p4scale
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = p5scale
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = d1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = d2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = d3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = d4scale
            end do
            do k = i, npole
               kk = ipole(k)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
c
c     calculate the error function damping terms
c
                  if (r2 .le. cut2) then
                     r = sqrt(r2)
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) / r
                     alsq2 = 2.0d0 * aewald**2
                     alsq2n = 0.0d0
                     if (aewald .gt. 0.0d0)
     &                  alsq2n = 1.0d0 / (sqrtpi*aewald)
                     exp2a = exp(-ralpha**2)
                     do j = 1, 3
                        bfac = dble(j+j-1)
                        alsq2n = alsq2 * alsq2n
                        bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
                     end do
c
c     compute the error function scaled and unscaled terms
c
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                           scale7 = 1.0d0 - expdamp
     &                                 *(1.0d0-damp+0.6d0*damp**2)
                        end if
                     end if
                     dsc3 = scale3
                     dsc5 = scale5
                     dsc7 = scale7
                     psc3 = scale3
                     psc5 = scale5
                     psc7 = scale7
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           dsc3 = scale3 * dscale(kk)
                           dsc5 = scale5 * dscale(kk)
                           dsc7 = scale7 * dscale(kk)
                           psc3 = scale3 * pscale(kk)
                           psc5 = scale5 * pscale(kk)
                           psc7 = scale7 * pscale(kk)
                        end if
                     end if
                     drr3 = (1.0d0-dsc3) / (r*r2)
                     drr5 = 3.0d0 * (1.0d0-dsc5) / (r*r2*r2)
                     drr7 = 15.0d0 * (1.0d0-dsc7) / (r*r2*r2*r2)
                     prr3 = (1.0d0-psc3) / (r*r2)
                     prr5 = 3.0d0 * (1.0d0-psc5) / (r*r2*r2)
                     prr7 = 15.0d0 * (1.0d0-psc7) / (r*r2*r2*r2)
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
                     fim(1) = -xr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                           - bn(1)*dkx + 2.0d0*bn(2)*qkx
                     fim(2) = -yr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                           - bn(1)*dky + 2.0d0*bn(2)*qky
                     fim(3) = -zr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                           - bn(1)*dkz + 2.0d0*bn(2)*qkz
                     fkm(1) = xr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                           - bn(1)*dix - 2.0d0*bn(2)*qix
                     fkm(2) = yr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                           - bn(1)*diy - 2.0d0*bn(2)*qiy
                     fkm(3) = zr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                           - bn(1)*diz - 2.0d0*bn(2)*qiz
                     fid(1) = -xr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                           - drr3*dkx + 2.0d0*drr5*qkx
                     fid(2) = -yr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                           - drr3*dky + 2.0d0*drr5*qky
                     fid(3) = -zr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                           - drr3*dkz + 2.0d0*drr5*qkz
                     fkd(1) = xr*(drr3*ci+drr5*dir+drr7*qir)
     &                           - drr3*dix - 2.0d0*drr5*qix
                     fkd(2) = yr*(drr3*ci+drr5*dir+drr7*qir)
     &                           - drr3*diy - 2.0d0*drr5*qiy
                     fkd(3) = zr*(drr3*ci+drr5*dir+drr7*qir)
     &                           - drr3*diz - 2.0d0*drr5*qiz
                     fip(1) = -xr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                           - prr3*dkx + 2.0d0*prr5*qkx
                     fip(2) = -yr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                           - prr3*dky + 2.0d0*prr5*qky
                     fip(3) = -zr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                           - prr3*dkz + 2.0d0*prr5*qkz
                     fkp(1) = xr*(prr3*ci+prr5*dir+prr7*qir)
     &                           - prr3*dix - 2.0d0*prr5*qix
                     fkp(2) = yr*(prr3*ci+prr5*dir+prr7*qir)
     &                           - prr3*diy - 2.0d0*prr5*qiy
                     fkp(3) = zr*(prr3*ci+prr5*dir+prr7*qir)
     &                           - prr3*diz - 2.0d0*prr5*qiz
c
c     increment the field at each site due to this interaction
c
                     do j = 1, 3
                        field(j,i) = field(j,i) + fim(j) - fid(j)
                        fieldp(j,i) = fieldp(j,i) + fim(j) - fip(j)
                        if (ii .ne. kk) then
                           field(j,k) = field(j,k) + fkm(j) - fkd(j)
                           fieldp(j,k) = fieldp(j,k) + fkm(j) - fkp(j)
                        end if
                     end do
                  end if
               end do
            end do
c
c     reset interaction scaling coefficients for connected atoms
c
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = 1.0d0
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = 1.0d0
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = 1.0d0
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = 1.0d0
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine udirect2b  --  Ewald real direct field via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "udirect2b" computes the real space contribution of the permanent
c     atomic multipole moments to the field via a neighbor list
c
c
      subroutine udirect2b (field,fieldp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'couple.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      integer i,j,k
      integer ii,kk,kkk
      real*8 xr,yr,zr,r,r2
      real*8 erfc,bfac,exp2a
      real*8 drr3,drr5,drr7
      real*8 prr3,prr5,prr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qiyy,qizz
      real*8 qixy,qixz,qiyz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkyy,qkzz
      real*8 qkxy,qkxz,qkyz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 ralpha
      real*8 alsq2,alsq2n
      real*8 pdi,pti,pgamma
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 bn(0:3)
      real*8 fim(3),fkm(3)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      character*6 mode
      external erfc
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     initialize local variables for OpenMP calculation
c
      do i = 1, npole
         do j = 1, 3
            fieldt(j,i) = 0.0d0
            fieldtp(j,i) = 0.0d0
         end do
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,pdi,pti,
!$OMP& ci,dix,diy,diz,qixx,qixy,qixz,qiyy,qiyz,qizz,kkk,kk,xr,yr,zr,
!$OMP& r2,r,ck,dkx,dky,dkz,qkxx,qkxy,qkxz,qkyy,qkyz,qkzz,ralpha,bn,
!$OMP& alsq2,alsq2n,exp2a,bfac,scale3,scale5,scale7,damp,pgamma,
!$OMP& dsc3,dsc5,dsc7,psc3,psc5,psc7,drr3,drr5,drr7,prr3,prr5,prr7,
!$OMP& dir,qix,qiy,qiz,qir,dkr,qkx,qky,qkz,qkr,fim,fkm,fid,fkd,fip,
!$OMP& fkp,expdamp) firstprivate(dscale,pscale)
!$OMP DO reduction(+:fieldt,fieldtp) schedule(dynamic)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 3
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     compute the error function scaled and unscaled terms
c
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                     scale7 = 1.0d0 - expdamp
     &                           *(1.0d0-damp+0.6d0*damp**2)
                  end if
               end if
               dsc3 = scale3 * dscale(kk)
               dsc5 = scale5 * dscale(kk)
               dsc7 = scale7 * dscale(kk)
               psc3 = scale3 * pscale(kk)
               psc5 = scale5 * pscale(kk)
               psc7 = scale7 * pscale(kk)
               drr3 = (1.0d0-dsc3) / (r*r2)
               drr5 = 3.0d0 * (1.0d0-dsc5) / (r*r2*r2)
               drr7 = 15.0d0 * (1.0d0-dsc7) / (r*r2*r2*r2)
               prr3 = (1.0d0-psc3) / (r*r2)
               prr5 = 3.0d0 * (1.0d0-psc5) / (r*r2*r2)
               prr7 = 15.0d0 * (1.0d0-psc7) / (r*r2*r2*r2)
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
               fim(1) = -xr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dkx + 2.0d0*bn(2)*qkx
               fim(2) = -yr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dky + 2.0d0*bn(2)*qky
               fim(3) = -zr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dkz + 2.0d0*bn(2)*qkz
               fkm(1) = xr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*dix - 2.0d0*bn(2)*qix
               fkm(2) = yr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*diy - 2.0d0*bn(2)*qiy
               fkm(3) = zr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*diz - 2.0d0*bn(2)*qiz
               fid(1) = -xr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dkx + 2.0d0*drr5*qkx
               fid(2) = -yr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dky + 2.0d0*drr5*qky
               fid(3) = -zr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dkz + 2.0d0*drr5*qkz
               fkd(1) = xr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*dix - 2.0d0*drr5*qix
               fkd(2) = yr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*diy - 2.0d0*drr5*qiy
               fkd(3) = zr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*diz - 2.0d0*drr5*qiz
               fip(1) = -xr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dkx + 2.0d0*prr5*qkx
               fip(2) = -yr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dky + 2.0d0*prr5*qky
               fip(3) = -zr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dkz + 2.0d0*prr5*qkz
               fkp(1) = xr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*dix - 2.0d0*prr5*qix
               fkp(2) = yr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*diy - 2.0d0*prr5*qiy
               fkp(3) = zr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*diz - 2.0d0*prr5*qiz
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  fieldt(j,i) = fieldt(j,i) + fim(j) - fid(j)
                  fieldt(j,k) = fieldt(j,k) + fkm(j) - fkd(j)
                  fieldtp(j,i) = fieldtp(j,i) + fim(j) - fip(j)
                  fieldtp(j,k) = fieldtp(j,k) + fkm(j) - fkp(j)
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
!$OMP END DO
c
c     end OpenMP directives for the major loop structure
c
!$OMP DO
      do i = 1, npole
         do j = 1, 3
            field(j,i) = fieldt(j,i) + field(j,i)
            fieldp(j,i) = fieldtp(j,i) + fieldp(j,i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine umutual1  --  Ewald recip mutual induced field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "umutual1" computes the reciprocal space contribution of the
c     induced atomic dipole moments to the field
c
c
      subroutine umutual1 (field,fieldp)
      implicit none
      include 'sizes.i'
      include 'boxes.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
      include 'pme.i'
      include 'polar.i'
      integer i,j,k
      real*8 term
      real*8 a(3,3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
      real*8, allocatable :: fdip_phi1(:,:)
      real*8, allocatable :: fdip_phi2(:,:)
      real*8, allocatable :: fdip_sum_phi(:,:)
      real*8, allocatable :: dipfield1(:,:)
      real*8, allocatable :: dipfield2(:,:)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (fuind(3,npole))
      allocate (fuinp(3,npole))
      allocate (fdip_phi1(10,npole))
      allocate (fdip_phi2(10,npole))
      allocate (fdip_sum_phi(20,npole))
      allocate (dipfield1(3,npole))
      allocate (dipfield2(3,npole))
c
c     convert Cartesian dipoles to fractional coordinates
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
      do i = 1, npole
         do k = 1, 3
            fuind(k,i) = a(k,1)*uind(1,i) + a(k,2)*uind(2,i)
     &                      + a(k,3)*uind(3,i)
            fuinp(k,i) = a(k,1)*uinp(1,i) + a(k,2)*uinp(2,i)
     &                      + a(k,3)*uinp(3,i)
         end do
      end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_uind (fuind,fuinp)
      call fftfront
c
c     complete the transformation of the PME grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get field
c
      call fftback
      call fphi_uind (fdip_phi1,fdip_phi2,fdip_sum_phi)
c
c     convert the dipole fields from fractional to Cartesian
c
      do i = 1, 3
         a(i,1) = dble(nfft1) * recip(i,1)
         a(i,2) = dble(nfft2) * recip(i,2)
         a(i,3) = dble(nfft3) * recip(i,3)
      end do
      do i = 1, npole
         do k = 1, 3
            dipfield1(k,i) = a(k,1)*fdip_phi1(2,i)
     &                          + a(k,2)*fdip_phi1(3,i)
     &                          + a(k,3)*fdip_phi1(4,i)
            dipfield2(k,i) = a(k,1)*fdip_phi2(2,i)
     &                          + a(k,2)*fdip_phi2(3,i)
     &                          + a(k,3)*fdip_phi2(4,i)
         end do
      end do
c
c     increment the field at each multipole site
c
      do i = 1, npole
         do k = 1, 3
            field(k,i) = field(k,i) - dipfield1(k,i)
            fieldp(k,i) = fieldp(k,i) - dipfield2(k,i)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (fdip_phi1)
      deallocate (fdip_phi2)
      deallocate (fdip_sum_phi)
      deallocate (dipfield1)
      deallocate (dipfield2)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine umutual2a  --  Ewald real mutual field via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "umutual2a" computes the real space contribution of the induced
c     atomic dipole moments to the field via a double loop
c
c
      subroutine umutual2a (field,fieldp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 erfc,bfac,exp2a
      real*8 duir,dukr
      real*8 puir,pukr
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 ralpha
      real*8 alsq2,alsq2n
      real*8 pdi,pti,pgamma
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 bn(0:2)
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      character*6 mode
      external erfc
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
c
c     set array needed to scale connected atom interactions
c
      do i = 1, n
         dscale(i) = 1.0d0
      end do
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         duix = uind(1,i)
         duiy = uind(2,i)
         duiz = uind(3,i)
         puix = uinp(1,i)
         puiy = uinp(2,i)
         puiz = uinp(3,i)
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = u4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               dukx = uind(1,k)
               duky = uind(2,k)
               dukz = uind(3,k)
               pukx = uinp(1,k)
               puky = uinp(2,k)
               pukz = uinp(3,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 2
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     compute the error function scaled and unscaled terms
c
               scale3 = dscale(kk)
               scale5 = dscale(kk)
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = scale3 * (1.0d0-expdamp)
                     scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
                  end if
               end if
               rr3 = (1.0d0-scale3) / (r*r2)
               rr5 = 3.0d0 * (1.0d0-scale5) / (r*r2*r2)
               duir = xr*duix + yr*duiy + zr*duiz
               dukr = xr*dukx + yr*duky + zr*dukz
               puir = xr*puix + yr*puiy + zr*puiz
               pukr = xr*pukx + yr*puky + zr*pukz
               fimd(1) = -bn(1)*dukx + bn(2)*dukr*xr
               fimd(2) = -bn(1)*duky + bn(2)*dukr*yr
               fimd(3) = -bn(1)*dukz + bn(2)*dukr*zr
               fkmd(1) = -bn(1)*duix + bn(2)*duir*xr
               fkmd(2) = -bn(1)*duiy + bn(2)*duir*yr
               fkmd(3) = -bn(1)*duiz + bn(2)*duir*zr
               fimp(1) = -bn(1)*pukx + bn(2)*pukr*xr
               fimp(2) = -bn(1)*puky + bn(2)*pukr*yr
               fimp(3) = -bn(1)*pukz + bn(2)*pukr*zr
               fkmp(1) = -bn(1)*puix + bn(2)*puir*xr
               fkmp(2) = -bn(1)*puiy + bn(2)*puir*yr
               fkmp(3) = -bn(1)*puiz + bn(2)*puir*zr
               fid(1) = -rr3*dukx + rr5*dukr*xr
               fid(2) = -rr3*duky + rr5*dukr*yr
               fid(3) = -rr3*dukz + rr5*dukr*zr
               fkd(1) = -rr3*duix + rr5*duir*xr
               fkd(2) = -rr3*duiy + rr5*duir*yr
               fkd(3) = -rr3*duiz + rr5*duir*zr
               fip(1) = -rr3*pukx + rr5*pukr*xr
               fip(2) = -rr3*puky + rr5*pukr*yr
               fip(3) = -rr3*pukz + rr5*pukr*zr
               fkp(1) = -rr3*puix + rr5*puir*xr
               fkp(2) = -rr3*puiy + rr5*puir*yr
               fkp(3) = -rr3*puiz + rr5*puir*zr
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  field(j,i) = field(j,i) + fimd(j) - fid(j)
                  field(j,k) = field(j,k) + fkmd(j) - fkd(j)
                  fieldp(j,i) = fieldp(j,i) + fimp(j) - fip(j)
                  fieldp(j,k) = fieldp(j,k) + fkmp(j) - fkp(j)
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do i = 1, npole
            ii = ipole(i)
            pdi = pdamp(i)
            pti = thole(i)
            duix = uind(1,i)
            duiy = uind(2,i)
            duiz = uind(3,i)
            puix = uinp(1,i)
            puiy = uinp(2,i)
            puiz = uinp(3,i)
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = u1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = u2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = u3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = u4scale
            end do
            do k = i, npole
               kk = ipole(k)
               dukx = uind(1,k)
               duky = uind(2,k)
               dukz = uind(3,k)
               pukx = uinp(1,k)
               puky = uinp(2,k)
               pukz = uinp(3,k)
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
c
c     calculate the error function damping terms
c
                  if (r2 .le. cut2) then
                     r = sqrt(r2)
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) / r
                     alsq2 = 2.0d0 * aewald**2
                     alsq2n = 0.0d0
                     if (aewald .gt. 0.0d0)
     &                  alsq2n = 1.0d0 / (sqrtpi*aewald)
                     exp2a = exp(-ralpha**2)
                     do j = 1, 2
                        bfac = dble(j+j-1)
                        alsq2n = alsq2 * alsq2n
                        bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
                     end do
c
c     compute the error function scaled and unscaled terms
c
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                        end if
                     end if
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           scale3 = scale3 * dscale(kk)
                           scale5 = scale5 * dscale(kk)
                        end if
                     end if
                     rr3 = (1.0d0-scale3) / (r*r2)
                     rr5 = 3.0d0 * (1.0d0-scale5) / (r*r2*r2)
                     duir = xr*duix + yr*duiy + zr*duiz
                     dukr = xr*dukx + yr*duky + zr*dukz
                     puir = xr*puix + yr*puiy + zr*puiz
                     pukr = xr*pukx + yr*puky + zr*pukz
                     fimd(1) = -bn(1)*dukx + bn(2)*dukr*xr
                     fimd(2) = -bn(1)*duky + bn(2)*dukr*yr
                     fimd(3) = -bn(1)*dukz + bn(2)*dukr*zr
                     fkmd(1) = -bn(1)*duix + bn(2)*duir*xr
                     fkmd(2) = -bn(1)*duiy + bn(2)*duir*yr
                     fkmd(3) = -bn(1)*duiz + bn(2)*duir*zr
                     fimp(1) = -bn(1)*pukx + bn(2)*pukr*xr
                     fimp(2) = -bn(1)*puky + bn(2)*pukr*yr
                     fimp(3) = -bn(1)*pukz + bn(2)*pukr*zr
                     fkmp(1) = -bn(1)*puix + bn(2)*puir*xr
                     fkmp(2) = -bn(1)*puiy + bn(2)*puir*yr
                     fkmp(3) = -bn(1)*puiz + bn(2)*puir*zr
                     fid(1) = -rr3*dukx + rr5*dukr*xr
                     fid(2) = -rr3*duky + rr5*dukr*yr
                     fid(3) = -rr3*dukz + rr5*dukr*zr
                     fkd(1) = -rr3*duix + rr5*duir*xr
                     fkd(2) = -rr3*duiy + rr5*duir*yr
                     fkd(3) = -rr3*duiz + rr5*duir*zr
                     fip(1) = -rr3*pukx + rr5*pukr*xr
                     fip(2) = -rr3*puky + rr5*pukr*yr
                     fip(3) = -rr3*pukz + rr5*pukr*zr
                     fkp(1) = -rr3*puix + rr5*puir*xr
                     fkp(2) = -rr3*puiy + rr5*puir*yr
                     fkp(3) = -rr3*puiz + rr5*puir*zr
c
c     increment the field at each site due to this interaction
c
                     do j = 1, 3
                        field(j,i) = field(j,i) + fimd(j) - fid(j)
                        fieldp(j,i) = fieldp(j,i) + fimp(j) - fip(j)
                        if (ii .ne. kk) then
                           field(j,k) = field(j,k) + fkmd(j) - fkd(j)
                           fieldp(j,k) = fieldp(j,k) + fkmp(j) - fkp(j)
                        end if
                     end do
                  end if
               end do
            end do
c
c     reset interaction scaling coefficients for connected atoms
c
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = 1.0d0
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = 1.0d0
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = 1.0d0
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine umutual2b  --  Ewald real mutual field via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "umutual2b" computes the real space contribution of the induced
c     atomic dipole moments to the field via a neighbor list
c
c
      subroutine umutual2b (field,fieldp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'couple.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      integer i,j,k
      integer ii,kk,kkk
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 erfc,bfac,exp2a
      real*8 duir,dukr
      real*8 puir,pukr
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 ralpha
      real*8 alsq2,alsq2n
      real*8 pdi,pti,pgamma
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 bn(0:2)
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      character*6 mode
      external erfc
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
c
c     set array needed to scale connected atom interactions
c
      do i = 1, n
         dscale(i) = 1.0d0
      end do
c
c     initialize local variables for OpenMP calculation
c
      do i = 1, npole
         do j = 1, 3
            fieldt(j,i) = 0.0d0
            fieldtp(j,i) = 0.0d0
         end do
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(xr,yr,zr,r,r2,rr3,rr5,
!$OMP& bfac,exp2a,duir,dukr,puir,pukr,pdi,pti,expdamp,
!$OMP& duix,duiy,duiz,puix,puiy,puiz,dukx,duky,dukz,pukx,puky,pukz,
!$OMP& ralpha,damp,alsq2,alsq2n,scale3,scale5,bn,fimd,fkmd,
!$OMP& fimp,fkmp,fid,fkd,fip,fkp,i,j,k,ii,kk,kkk)
!$OMP& firstprivate(dscale)
!$OMP DO reduction(+:fieldt,fieldtp) schedule(dynamic)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         duix = uind(1,i)
         duiy = uind(2,i)
         duiz = uind(3,i)
         puix = uinp(1,i)
         puiy = uinp(2,i)
         puiz = uinp(3,i)
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = u4scale
         end do
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               dukx = uind(1,k)
               duky = uind(2,k)
               dukz = uind(3,k)
               pukx = uinp(1,k)
               puky = uinp(2,k)
               pukz = uinp(3,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 2
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     compute the error function scaled and unscaled terms
c
               scale3 = dscale(kk)
               scale5 = dscale(kk)
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = scale3 * (1.0d0-expdamp)
                     scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
                  end if
               end if
               rr3 = (1.0d0-scale3) / (r*r2)
               rr5 = 3.0d0 * (1.0d0-scale5) / (r*r2*r2)
               duir = xr*duix + yr*duiy + zr*duiz
               dukr = xr*dukx + yr*duky + zr*dukz
               puir = xr*puix + yr*puiy + zr*puiz
               pukr = xr*pukx + yr*puky + zr*pukz
               fimd(1) = -bn(1)*dukx + bn(2)*dukr*xr
               fimd(2) = -bn(1)*duky + bn(2)*dukr*yr
               fimd(3) = -bn(1)*dukz + bn(2)*dukr*zr
               fkmd(1) = -bn(1)*duix + bn(2)*duir*xr
               fkmd(2) = -bn(1)*duiy + bn(2)*duir*yr
               fkmd(3) = -bn(1)*duiz + bn(2)*duir*zr
               fimp(1) = -bn(1)*pukx + bn(2)*pukr*xr
               fimp(2) = -bn(1)*puky + bn(2)*pukr*yr
               fimp(3) = -bn(1)*pukz + bn(2)*pukr*zr
               fkmp(1) = -bn(1)*puix + bn(2)*puir*xr
               fkmp(2) = -bn(1)*puiy + bn(2)*puir*yr
               fkmp(3) = -bn(1)*puiz + bn(2)*puir*zr
               fid(1) = -rr3*dukx + rr5*dukr*xr
               fid(2) = -rr3*duky + rr5*dukr*yr
               fid(3) = -rr3*dukz + rr5*dukr*zr
               fkd(1) = -rr3*duix + rr5*duir*xr
               fkd(2) = -rr3*duiy + rr5*duir*yr
               fkd(3) = -rr3*duiz + rr5*duir*zr
               fip(1) = -rr3*pukx + rr5*pukr*xr
               fip(2) = -rr3*puky + rr5*pukr*yr
               fip(3) = -rr3*pukz + rr5*pukr*zr
               fkp(1) = -rr3*puix + rr5*puir*xr
               fkp(2) = -rr3*puiy + rr5*puir*yr
               fkp(3) = -rr3*puiz + rr5*puir*zr
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  fieldt(j,i) = fieldt(j,i) + fimd(j) - fid(j)
                  fieldt(j,k) = fieldt(j,k) + fkmd(j) - fkd(j)
                  fieldtp(j,i) = fieldtp(j,i) + fimp(j) - fip(j)
                  fieldtp(j,k) = fieldtp(j,k) + fkmp(j) - fkp(j)
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
!$OMP END DO
c
c     end OpenMP directives for the major loop structure
c
!$OMP DO
      do i = 1, npole
         do j = 1, 3
            field(j,i) = fieldt(j,i) + field(j,i)
            fieldp(j,i) = fieldtp(j,i) + fieldp(j,i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine induce0e  --  Kirkwood SCRF induced dipoles  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "induce0e" computes the induced dipole moments at polarizable
c     sites for generalized Kirkwood SCRF and vacuum environments
c
c
      subroutine induce0e
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'couple.i'
      include 'gkstuf.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'solute.i'
      include 'units.i'
      integer i,j,k
      integer ii,kk
      integer iter,maxiter
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,uxi,uyi,uzi
      real*8 qxxi,qxyi,qxzi
      real*8 qyyi,qyzi,qzzi
      real*8 ck,uxk,uyk,uzk
      real*8 qxxk,qxyk,qxzk
      real*8 qyyk,qyzk,qzzk
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 dir,duir,puir
      real*8 dkr,dukr,pukr
      real*8 duixs,duiys,duizs
      real*8 puixs,puiys,puizs
      real*8 dukxs,dukys,dukzs
      real*8 pukxs,pukys,pukzs
      real*8 duirs,puirs
      real*8 dukrs,pukrs
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 epsds,epsps
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 rb2,rbi,rbk
      real*8 dwater,fc,fd,fq
      real*8 gf,gf2,gf3,gf5,gf7
      real*8 expterm,expc,expc1
      real*8 dexpc,expcdexpc
      real*8 a(0:3,0:2)
      real*8 gc(4),gux(10)
      real*8 guy(10),guz(10)
      real*8 gqxx(4),gqxy(4)
      real*8 gqxz(4),gqyy(4)
      real*8 gqyz(4),gqzz(4)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 fids(3),fkds(3)
      real*8 fips(3),fkps(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: fields(:,:)
      real*8, allocatable :: fieldps(:,:)
      real*8, allocatable :: gkfield(:,:)
      real*8, allocatable :: gkfieldp(:,:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: udirs(:,:)
      real*8, allocatable :: udirps(:,:)
      real*8, allocatable :: uold(:,:)
      real*8, allocatable :: uoldp(:,:)
      real*8, allocatable :: uolds(:,:)
      real*8, allocatable :: uoldps(:,:)
      logical proceed,done
      character*6 mode
c
c
c     zero out the induced dipoles at each site; uind & uinp are
c     vacuum dipoles, uinds & uinps are SCRF dipoles
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
            uinds(j,i) = 0.0d0
            uinps(j,i) = 0.0d0
         end do
      end do
      if (.not.use_polar .and. .not.use_solv)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
      allocate (gkfield(3,npole))
      allocate (gkfieldp(3,npole))
c
c     zero out the fields at each site; field & fieldp are solute
c     fields, gkfield & gkfieldp are reaction fields
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
            gkfield(j,i) = 0.0d0
            gkfieldp(j,i) = 0.0d0
         end do
      end do
      dwater = 78.3d0
      fc = 1.0d0 * (1.0d0 - dwater)/(0.0d0 + 1.0d0*dwater)
      fd = 2.0d0 * (1.0d0 - dwater)/(1.0d0 + 2.0d0*dwater)
      fq = 3.0d0 * (1.0d0 - dwater)/(2.0d0 + 3.0d0*dwater)
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     compute the direct induced dipole moment at each atom, and
c     another set that also includes RF due to permanent multipoles
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         rbi = rborn(ii)
         ci = rpole(1,i)
         uxi = rpole(2,i)
         uyi = rpole(3,i)
         uzi = rpole(4,i)
         qxxi = rpole(5,i)
         qxyi = rpole(6,i)
         qxzi = rpole(7,i)
         qyyi = rpole(9,i)
         qyzi = rpole(10,i)
         qzzi = rpole(13,i)
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do k = i, npole
            kk = ipole(k)
            rbk = rborn(kk)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  uxk = rpole(2,k)
                  uyk = rpole(3,k)
                  uzk = rpole(4,k)
                  qxxk = rpole(5,k)
                  qxyk = rpole(6,k)
                  qxzk = rpole(7,k)
                  qyyk = rpole(9,k)
                  qyzk = rpole(10,k)
                  qzzk = rpole(13,k)
c
c     self-interactions for the solute field are skipped
c
                  if (i .ne. k) then
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                           scale7 = 1.0d0 - expdamp
     &                                 *(1.0d0-damp+0.6d0*damp**2)
                        end if
                     end if
                     rr3 = scale3 / (r*r2)
                     rr5 = 3.0d0 * scale5 / (r*r2*r2)
                     rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                     dir = uxi*xr + uyi*yr + uzi*zr
                     qix = qxxi*xr + qxyi*yr + qxzi*zr
                     qiy = qxyi*xr + qyyi*yr + qyzi*zr
                     qiz = qxzi*xr + qyzi*yr + qzzi*zr
                     qir = qix*xr + qiy*yr + qiz*zr
                     dkr = uxk*xr + uyk*yr + uzk*zr
                     qkx = qxxk*xr + qxyk*yr + qxzk*zr
                     qky = qxyk*xr + qyyk*yr + qyzk*zr
                     qkz = qxzk*xr + qyzk*yr + qzzk*zr
                     qkr = qkx*xr + qky*yr + qkz*zr
                     fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*uxk + 2.0d0*rr5*qkx
                     fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*uyk + 2.0d0*rr5*qky
                     fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*uzk + 2.0d0*rr5*qkz
                     fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*uxi - 2.0d0*rr5*qix
                     fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*uyi - 2.0d0*rr5*qiy
                     fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*uzi - 2.0d0*rr5*qiz
                     do j = 1, 3
                        field(j,i) = field(j,i) + fid(j)*dscale(kk)
                        field(j,k) = field(j,k) + fkd(j)*dscale(kk)
                        fieldp(j,i) = fieldp(j,i) + fid(j)*pscale(kk)
                        fieldp(j,k) = fieldp(j,k) + fkd(j)*pscale(kk)
                     end do
                  end if
                  rb2 = rbi * rbk
                  expterm = exp(-r2/(gkc*rb2))
                  expc = expterm / gkc
                  dexpc = -2.0d0 / (gkc*rb2)
                  gf2 = 1.0d0 / (r2+rb2*expterm)
                  gf = sqrt(gf2)
                  gf3 = gf2 * gf
                  gf5 = gf3 * gf2
                  gf7 = gf5 * gf2
c
c     reaction potential auxiliary terms
c
                  a(0,0) = gf
                  a(1,0) = -gf3
                  a(2,0) = 3.0d0 * gf5
                  a(3,0) = -15.0d0 * gf7
c
c     reaction potential gradient auxiliary terms
c
                  expc1 = 1.0d0 - expc
                  a(0,1) = expc1 * a(1,0)
                  a(1,1) = expc1 * a(2,0)
                  a(2,1) = expc1 * a(3,0)
c
c     dipole second reaction potential gradient auxiliary term
c
                  expcdexpc = -expc * dexpc
                  a(1,2) = expc1*a(2,1) + expcdexpc*a(2,0)
c
c     multiply the auxillary terms by dielectric functions
c
                  a(0,1) = fc * a(0,1)
                  a(1,0) = fd * a(1,0)
                  a(1,1) = fd * a(1,1)
                  a(1,2) = fd * a(1,2)
                  a(2,0) = fq * a(2,0)
                  a(2,1) = fq * a(2,1)
c
c     unweighted dipole reaction potential tensor
c
                  gux(1) = xr * a(1,0)
                  guy(1) = yr * a(1,0)
                  guz(1) = zr * a(1,0)
c
c     unweighted reaction potential gradient tensor
c
                  gc(2) = xr * a(0,1)
                  gc(3) = yr * a(0,1)
                  gc(4) = zr * a(0,1)
                  gux(2) = a(1,0) + xr2*a(1,1)
                  gux(3) = xr * yr * a(1,1)
                  gux(4) = xr * zr * a(1,1)
                  guy(2) = gux(3)
                  guy(3) = a(1,0) + yr2*a(1,1)
                  guy(4) = yr * zr * a(1,1)
                  guz(2) = gux(4)
                  guz(3) = guy(4)
                  guz(4) = a(1,0) + zr2*a(1,1)
                  gqxx(2) = xr * (2.0d0*a(2,0)+xr2*a(2,1))
                  gqxx(3) = yr * xr2*a(2,1)
                  gqxx(4) = zr * xr2*a(2,1)
                  gqyy(2) = xr * yr2*a(2,1)
                  gqyy(3) = yr * (2.0d0*a(2,0)+yr2*a(2,1))
                  gqyy(4) = zr * yr2 * a(2,1)
                  gqzz(2) = xr * zr2 * a(2,1)
                  gqzz(3) = yr * zr2 * a(2,1)
                  gqzz(4) = zr * (2.0d0*a(2,0)+zr2*a(2,1))
                  gqxy(2) = yr * (a(2,0)+xr2*a(2,1))
                  gqxy(3) = xr * (a(2,0)+yr2*a(2,1))
                  gqxy(4) = zr * xr * yr * a(2,1)
                  gqxz(2) = zr * (a(2,0)+xr2*a(2,1))
                  gqxz(3) = gqxy(4)
                  gqxz(4) = xr * (a(2,0)+zr2*a(2,1))
                  gqyz(2) = gqxy(4)
                  gqyz(3) = zr * (a(2,0)+yr2*a(2,1))
                  gqyz(4) = yr * (a(2,0)+zr2*a(2,1))
c
c     unweighted dipole second reaction potential gradient tensor
c
                  gux(5) = xr * (3.0d0*a(1,1)+xr2*a(1,2))
                  gux(6) = yr * (a(1,1)+xr2*a(1,2))
                  gux(7) = zr * (a(1,1)+xr2*a(1,2))
                  gux(8) = xr * (a(1,1)+yr2*a(1,2))
                  gux(9) = zr * xr * yr * a(1,2)
                  gux(10) = xr * (a(1,1)+zr2*a(1,2))
                  guy(5) = yr * (a(1,1)+xr2*a(1,2))
                  guy(6) = xr * (a(1,1)+yr2*a(1,2))
                  guy(7) = gux(9)
                  guy(8) = yr * (3.0d0*a(1,1)+yr2*a(1,2))
                  guy(9) = zr * (a(1,1)+yr2*a(1,2))
                  guy(10) = yr * (a(1,1)+zr2*a(1,2))
                  guz(5) = zr * (a(1,1)+xr2*a(1,2))
                  guz(6) = gux(9)
                  guz(7) = xr * (a(1,1)+zr2*a(1,2))
                  guz(8) = zr * (a(1,1)+yr2*a(1,2))
                  guz(9) = yr * (a(1,1)+zr2*a(1,2))
                  guz(10) = zr * (3.0d0*a(1,1)+zr2*a(1,2))
c
c     generalized Kirkwood permanent reaction field
c
                  fid(1) = uxk*gux(2) + uyk*gux(3) + uzk*gux(4)
     &                        + 0.5d0 * (ck*gux(1) + qxxk*gux(5)
     &                            + qyyk*gux(8) + qzzk*gux(10)
     &                            + 2.0d0*(qxyk*gux(6)+qxzk*gux(7)
     &                                         +qyzk*gux(9)))
     &                        + 0.5d0 * (ck*gc(2) + qxxk*gqxx(2)
     &                            + qyyk*gqyy(2) + qzzk*gqzz(2)
     &                            + 2.0d0*(qxyk*gqxy(2)+qxzk*gqxz(2)
     &                                         +qyzk*gqyz(2)))
                  fid(2) = uxk*guy(2) + uyk*guy(3) + uzk*guy(4)
     &                        + 0.5d0 * (ck*guy(1) + qxxk*guy(5)
     &                            + qyyk*guy(8) + qzzk*guy(10)
     &                            + 2.0d0*(qxyk*guy(6)+qxzk*guy(7)
     &                                         +qyzk*guy(9)))
     &                        + 0.5d0 * (ck*gc(3) + qxxk*gqxx(3)
     &                            + qyyk*gqyy(3) + qzzk*gqzz(3)
     &                            + 2.0d0*(qxyk*gqxy(3)+qxzk*gqxz(3)
     &                                         +qyzk*gqyz(3)))
                  fid(3) = uxk*guz(2) + uyk*guz(3) + uzk*guz(4)
     &                        + 0.5d0 * (ck*guz(1) + qxxk*guz(5)
     &                            + qyyk*guz(8) + qzzk*guz(10)
     &                            + 2.0d0*(qxyk*guz(6)+qxzk*guz(7)
     &                                         +qyzk*guz(9)))
     &                        + 0.5d0 * (ck*gc(4) + qxxk*gqxx(4)
     &                            + qyyk*gqyy(4) + qzzk*gqzz(4)
     &                            + 2.0d0*(qxyk*gqxy(4)+qxzk*gqxz(4)
     &                                         +qyzk*gqyz(4)))
                  fkd(1) = uxi*gux(2) + uyi*gux(3) + uzi*gux(4)
     &                        - 0.5d0 * (ci*gux(1) + qxxi*gux(5)
     &                            + qyyi*gux(8) + qzzi*gux(10)
     &                            + 2.0d0*(qxyi*gux(6)+qxzi*gux(7)
     &                                         +qyzi*gux(9)))
     &                        - 0.5d0 * (ci*gc(2) + qxxi*gqxx(2)
     &                            + qyyi*gqyy(2) + qzzi*gqzz(2)
     &                            + 2.0d0*(qxyi*gqxy(2)+qxzi*gqxz(2)
     &                                         +qyzi*gqyz(2)))
                  fkd(2) = uxi*guy(2) + uyi*guy(3) + uzi*guy(4)
     &                        - 0.5d0 * (ci*guy(1) + qxxi*guy(5)
     &                            + qyyi*guy(8) + qzzi*guy(10)
     &                            + 2.0d0*(qxyi*guy(6)+qxzi*guy(7)
     &                                         +qyzi*guy(9)))
     &                        - 0.5d0 * (ci*gc(3) + qxxi*gqxx(3)
     &                            + qyyi*gqyy(3) + qzzi*gqzz(3)
     &                            + 2.0d0*(qxyi*gqxy(3)+qxzi*gqxz(3)
     &                                         +qyzi*gqyz(3)))
                  fkd(3) = uxi*guz(2) + uyi*guz(3) + uzi*guz(4)
     &                        - 0.5d0 * (ci*guz(1) + qxxi*guz(5)
     &                            + qyyi*guz(8) + qzzi*guz(10)
     &                            + 2.0d0*(qxyi*guz(6)+qxzi*guz(7)
     &                                         +qyzi*guz(9)))
     &                        - 0.5d0 * (ci*gc(4) + qxxi*gqxx(4)
     &                            + qyyi*gqyy(4) + qzzi*gqzz(4)
     &                            + 2.0d0*(qxyi*gqxy(4)+qxzi*gqxz(4)
     &                                         +qyzi*gqyz(4)))
c
c     scale the self-field by half, such that it sums to one below
c
                  if (i .eq. k) then
                     do j = 1,3
                        fid(j) = 0.5d0 * fid(j)
                        fkd(j) = 0.5d0 * fkd(j)
                     end do
                  end if
                  do j = 1, 3
                     gkfield(j,i) = gkfield(j,i) + fid(j)
                     gkfield(j,k) = gkfield(j,k) + fkd(j)
                     gkfieldp(j,i) = gkfieldp(j,i) + fid(j)
                     gkfieldp(j,k) = gkfieldp(j,k) + fkd(j)
                  end do
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (fields(3,npole))
      allocate (fieldps(3,npole))
      allocate (udir(3,npole))
      allocate (udirp(3,npole))
      allocate (udirs(3,npole))
      allocate (udirps(3,npole))
      allocate (uold(3,npole))
      allocate (uoldp(3,npole))
      allocate (uolds(3,npole))
      allocate (uoldps(3,npole))
c
c     set vacuum induced dipoles to polarizability times direct field;
c     set SCRF induced dipoles to polarizability times direct field
c     plus the GK reaction field due to permanent multipoles
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            udirs(j,i) = polarity(i) * (field(j,i)+gkfield(j,i))
            udirps(j,i) = polarity(i) * (fieldp(j,i)+gkfieldp(j,i))
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
            uinds(j,i) = udirs(j,i)
            uinps(j,i) = udirps(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         eps = 100.0d0
c
c     compute mutual induced dipole moments by an iterative method
c
         do while (.not. done)
            do i = 1, npole
               do j = 1, 3
                  field(j,i) = 0.0d0
                  fieldp(j,i) = 0.0d0
                  fields(j,i) = 0.0d0
                  fieldps(j,i) = 0.0d0
                  gkfield(j,i) = 0.0d0
                  gkfieldp(j,i) = 0.0d0
               end do
            end do
            do i = 1, npole
               ii = ipole(i)
               pdi = pdamp(i)
               pti = thole(i)
               rbi = rborn(ii)
               duix = uind(1,i)
               duiy = uind(2,i)
               duiz = uind(3,i)
               puix = uinp(1,i)
               puiy = uinp(2,i)
               puiz = uinp(3,i)
               duixs = uinds(1,i)
               duiys = uinds(2,i)
               duizs = uinds(3,i)
               puixs = uinps(1,i)
               puiys = uinps(2,i)
               puizs = uinps(3,i)
               do j = 1, np11(ii)
                  dscale(ip11(j,ii)) = u1scale
               end do
               do j = 1, np12(ii)
                  dscale(ip12(j,ii)) = u2scale
               end do
               do j = 1, np13(ii)
                  dscale(ip13(j,ii)) = u3scale
               end do
               do j = 1, np14(ii)
                  dscale(ip14(j,ii)) = u4scale
               end do
               do k = i, npole
                  kk = ipole(k)
                  rbk = rborn(kk)
                  proceed = .true.
                  if (use_intra)
     &               call groups (proceed,fgrp,ii,kk,0,0,0,0)
                  if (proceed) then
                     xr = x(kk) - x(ii)
                     yr = y(kk) - y(ii)
                     zr = z(kk) - z(ii)
                     xr2 = xr * xr
                     yr2 = yr * yr
                     zr2 = zr * zr
                     r2 = xr2 + yr2 + zr2
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        dukx = uind(1,k)
                        duky = uind(2,k)
                        dukz = uind(3,k)
                        pukx = uinp(1,k)
                        puky = uinp(2,k)
                        pukz = uinp(3,k)
                        dukxs = uinds(1,k)
                        dukys = uinds(2,k)
                        dukzs = uinds(3,k)
                        pukxs = uinps(1,k)
                        pukys = uinps(2,k)
                        pukzs = uinps(3,k)
                        if (i .ne. k) then
                           scale3 = dscale(kk)
                           scale5 = dscale(kk)
                           damp = pdi * pdamp(k)
                           if (damp .ne. 0.0d0) then
                              pgamma = min(pti,thole(k))
                              damp = -pgamma * (r/damp)**3
                              if (damp .gt. -50.0d0) then
                                 expdamp = exp(damp)
                                 scale3 = scale3 * (1.0d0-expdamp)
                                 scale5 = scale5 * (1.0d0-(1.0d0-damp)
     &                                                 *expdamp)
                              end if
                           end if
                           rr3 = scale3 / (r*r2)
                           rr5 = 3.0d0 * scale5 / (r*r2*r2)
                           duir = xr*duix + yr*duiy + zr*duiz
                           dukr = xr*dukx + yr*duky + zr*dukz
                           puir = xr*puix + yr*puiy + zr*puiz
                           pukr = xr*pukx + yr*puky + zr*pukz
                           duirs = xr*duixs + yr*duiys + zr*duizs
                           dukrs = xr*dukxs + yr*dukys + zr*dukzs
                           puirs = xr*puixs + yr*puiys + zr*puizs
                           pukrs = xr*pukxs + yr*pukys + zr*pukzs
                           fid(1) = -rr3*dukx + rr5*dukr*xr
                           fid(2) = -rr3*duky + rr5*dukr*yr
                           fid(3) = -rr3*dukz + rr5*dukr*zr
                           fkd(1) = -rr3*duix + rr5*duir*xr
                           fkd(2) = -rr3*duiy + rr5*duir*yr
                           fkd(3) = -rr3*duiz + rr5*duir*zr
                           fip(1) = -rr3*pukx + rr5*pukr*xr
                           fip(2) = -rr3*puky + rr5*pukr*yr
                           fip(3) = -rr3*pukz + rr5*pukr*zr
                           fkp(1) = -rr3*puix + rr5*puir*xr
                           fkp(2) = -rr3*puiy + rr5*puir*yr
                           fkp(3) = -rr3*puiz + rr5*puir*zr
                           fids(1) = -rr3*dukxs + rr5*dukrs*xr
                           fids(2) = -rr3*dukys + rr5*dukrs*yr
                           fids(3) = -rr3*dukzs + rr5*dukrs*zr
                           fkds(1) = -rr3*duixs + rr5*duirs*xr
                           fkds(2) = -rr3*duiys + rr5*duirs*yr
                           fkds(3) = -rr3*duizs + rr5*duirs*zr
                           fips(1) = -rr3*pukxs + rr5*pukrs*xr
                           fips(2) = -rr3*pukys + rr5*pukrs*yr
                           fips(3) = -rr3*pukzs + rr5*pukrs*zr
                           fkps(1) = -rr3*puixs + rr5*puirs*xr
                           fkps(2) = -rr3*puiys + rr5*puirs*yr
                           fkps(3) = -rr3*puizs + rr5*puirs*zr
                           do j = 1, 3
                              field(j,i) = field(j,i) + fid(j)
                              field(j,k) = field(j,k) + fkd(j)
                              fieldp(j,i) = fieldp(j,i) + fip(j)
                              fieldp(j,k) = fieldp(j,k) + fkp(j)
                              fields(j,i) = fields(j,i) + fids(j)
                              fields(j,k) = fields(j,k) + fkds(j)
                              fieldps(j,i) = fieldps(j,i) + fips(j)
                              fieldps(j,k) = fieldps(j,k) + fkps(j)
                           end do
                        end if
                        rb2 = rbi * rbk
                        expterm = exp(-r2/(gkc*rb2))
                        expc = expterm / gkc
                        dexpc = -2.0d0 / (gkc*rbi*rbk)
                        gf2 = 1.0d0 / (r2+rb2*expterm)
                        gf = sqrt(gf2)
                        gf3 = gf2 * gf
                        gf5 = gf3 * gf2
c
c     reaction potential auxiliary terms
c
                        a(1,0) = -gf3
                        a(2,0) = 3.0d0 * gf5
c
c     reaction potential gradient auxiliary terms
c
                        expc1 = 1.0d0 - expc
                        a(1,1) = expc1 * a(2,0)
c
c     unweighted dipole reaction potential gradient tensor
c
                        gux(2) = fd * (a(1,0) + xr2*a(1,1))
                        gux(3) = fd * xr*yr*a(1,1)
                        gux(4) = fd * xr*zr*a(1,1)
                        guy(2) = gux(3)
                        guy(3) = fd * (a(1,0) + yr2*a(1,1))
                        guy(4) = fd * yr*zr*a(1,1)
                        guz(2) = gux(4)
                        guz(3) = guy(4)
                        guz(4) = fd * (a(1,0) + zr2*a(1,1))
                        fids(1) = dukxs*gux(2)+dukys*guy(2)+dukzs*guz(2)
                        fids(2) = dukxs*gux(3)+dukys*guy(3)+dukzs*guz(3)
                        fids(3) = dukxs*gux(4)+dukys*guy(4)+dukzs*guz(4)
                        fkds(1) = duixs*gux(2)+duiys*guy(2)+duizs*guz(2)
                        fkds(2) = duixs*gux(3)+duiys*guy(3)+duizs*guz(3)
                        fkds(3) = duixs*gux(4)+duiys*guy(4)+duizs*guz(4)
                        fips(1) = pukxs*gux(2)+pukys*guy(2)+pukzs*guz(2)
                        fips(2) = pukxs*gux(3)+pukys*guy(3)+pukzs*guz(3)
                        fips(3) = pukxs*gux(4)+pukys*guy(4)+pukzs*guz(4)
                        fkps(1) = puixs*gux(2)+puiys*guy(2)+puizs*guz(2)
                        fkps(2) = puixs*gux(3)+puiys*guy(3)+puizs*guz(3)
                        fkps(3) = puixs*gux(4)+puiys*guy(4)+puizs*guz(4)
                        if (i .eq. k) then
                           do j = 1, 3
                              fids(j) = 0.5d0 * fids(j)
                              fkds(j) = 0.5d0 * fkds(j)
                              fips(j) = 0.5d0 * fips(j)
                              fkps(j) = 0.5d0 * fkps(j)
                           end do
                        end if
                        do j = 1, 3
                           gkfield(j,i) = gkfield(j,i) + fids(j)
                           gkfield(j,k) = gkfield(j,k) + fkds(j)
                           gkfieldp(j,i) = gkfieldp(j,i) + fips(j)
                           gkfieldp(j,k) = gkfieldp(j,k) + fkps(j)
                        end do
                     end if
                  end if
               end do
c
c     reset interaction scaling coefficients for connected atoms
c
               do j = 1, np11(ii)
                  dscale(ip11(j,ii)) = 1.0d0
               end do
               do j = 1, np12(ii)
                  dscale(ip12(j,ii)) = 1.0d0
               end do
               do j = 1, np13(ii)
                  dscale(ip13(j,ii)) = 1.0d0
               end do
               do j = 1, np14(ii)
                  dscale(ip14(j,ii)) = 1.0d0
               end do
            end do
c
c     check to see if the mutual induced dipoles have converged
c
            iter = iter + 1
            epsold = eps
            epsd = 0.0d0
            epsp = 0.0d0
            epsds = 0.0d0
            epsps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  uold(j,i) = uind(j,i)
                  uoldp(j,i) = uinp(j,i)
                  uolds(j,i) = uinds(j,i)
                  uoldps(j,i) = uinps(j,i)
                  uind(j,i) = udir(j,i) + polarity(i)*field(j,i)
                  uinp(j,i) = udirp(j,i) + polarity(i)*fieldp(j,i)
                  uinds(j,i) = udirs(j,i) + polarity(i)
     &                               *(fields(j,i)+gkfield(j,i))
                  uinps(j,i) = udirps(j,i) + polarity(i)
     &                                *(fieldps(j,i)+gkfieldp(j,i))
                  uind(j,i) = uold(j,i) + polsor*(uind(j,i)-uold(j,i))
                  uinp(j,i) = uoldp(j,i) + polsor*(uinp(j,i)-uoldp(j,i))
                  uinds(j,i) = uolds(j,i) + polsor
     &                              *(uinds(j,i)-uolds(j,i))
                  uinps(j,i) = uoldps(j,i) + polsor
     &                               *(uinps(j,i)-uoldps(j,i))
                  epsd = epsd + (uind(j,i)-uold(j,i))**2
                  epsp = epsp + (uinp(j,i)-uoldp(j,i))**2
                  epsds = epsds + (uinds(j,i)-uolds(j,i))**2
                  epsps = epsps + (uinps(j,i)-uoldps(j,i))**2
               end do
            end do
            eps = max(epsd,epsp,epsds,epsps)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. maxiter)  done = .true.
         end do
c
c     terminate the calculation if dipoles failed to converge
c
         if (eps .gt. poleps) then
            write (iout,30)
   30       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      deallocate (field)
      deallocate (fieldp)
      deallocate (fields)
      deallocate (fieldps)
      deallocate (gkfield)
      deallocate (gkfieldp)
      deallocate (udir)
      deallocate (udirp)
      deallocate (udirs)
      deallocate (udirps)
      deallocate (uold)
      deallocate (uoldp)
      deallocate (uolds)
      deallocate (uoldps)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine induce0f  --  Poisson-Boltzmann induced dipoles  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "induce0f" computes the induced dipole moments at polarizable
c     sites for Poisson-Boltzmann SCRF and vacuum environments
c
c
      subroutine induce0f
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'couple.i'
      include 'gkstuf.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'pbstuf.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'solute.i'
      include 'units.i'
      integer i,j,k
      integer ii,kk
      integer iter,maxiter
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,uxi,uyi,uzi
      real*8 qxxi,qxyi,qxzi
      real*8 qyyi,qyzi,qzzi
      real*8 ck,uxk,uyk,uzk
      real*8 qxxk,qxyk,qxzk
      real*8 qyyk,qyzk,qzzk
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 dir,duir,puir
      real*8 dkr,dukr,pukr
      real*8 duixs,duiys,duizs
      real*8 puixs,puiys,puizs
      real*8 dukxs,dukys,dukzs
      real*8 pukxs,pukys,pukzs
      real*8 duirs,puirs
      real*8 dukrs,pukrs
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 epsds,epsps
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 fids(3),fkds(3)
      real*8 fips(3),fkps(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: fields(:,:)
      real*8, allocatable :: fieldps(:,:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: udirs(:,:)
      real*8, allocatable :: udirps(:,:)
      real*8, allocatable :: uold(:,:)
      real*8, allocatable :: uoldp(:,:)
      real*8, allocatable :: uolds(:,:)
      real*8, allocatable :: uoldps(:,:)
      real*8, allocatable :: indpole(:,:)
      real*8, allocatable :: inppole(:,:)
      logical proceed,done
      character*6 mode
c
c
c     zero out the induced dipoles; uind and uinp are vacuum dipoles,
c     uinds and uinps are SCRF dipoles
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
            uinds(j,i) = 0.0d0
            uinps(j,i) = 0.0d0
         end do
      end do
      if (.not.use_polar .or. .not.use_solv)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
c
c     zero out the fields at each site; field and fieldp are
c     the solute fields
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     compute the direct induced dipole moment at each atom, and
c     another set that also includes RF due to permanent multipoles
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         uxi = rpole(2,i)
         uyi = rpole(3,i)
         uzi = rpole(4,i)
         qxxi = rpole(5,i)
         qxyi = rpole(6,i)
         qxzi = rpole(7,i)
         qyyi = rpole(9,i)
         qyzi = rpole(10,i)
         qzzi = rpole(13,i)
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  uxk = rpole(2,k)
                  uyk = rpole(3,k)
                  uzk = rpole(4,k)
                  qxxk = rpole(5,k)
                  qxyk = rpole(6,k)
                  qxzk = rpole(7,k)
                  qyyk = rpole(9,k)
                  qyzk = rpole(10,k)
                  qzzk = rpole(13,k)
c
c     self-interactions for the solute field are skipped
c
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                        scale7 = 1.0d0 - expdamp
     &                              *(1.0d0-damp+0.6d0*damp**2)
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                  dir = uxi*xr + uyi*yr + uzi*zr
                  qix = qxxi*xr + qxyi*yr + qxzi*zr
                  qiy = qxyi*xr + qyyi*yr + qyzi*zr
                  qiz = qxzi*xr + qyzi*yr + qzzi*zr
                  qir = qix*xr + qiy*yr + qiz*zr
                  dkr = uxk*xr + uyk*yr + uzk*zr
                  qkx = qxxk*xr + qxyk*yr + qxzk*zr
                  qky = qxyk*xr + qyyk*yr + qyzk*zr
                  qkz = qxzk*xr + qyzk*yr + qzzk*zr
                  qkr = qkx*xr + qky*yr + qkz*zr
                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*uxk + 2.0d0*rr5*qkx
                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*uyk + 2.0d0*rr5*qky
                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*uzk + 2.0d0*rr5*qkz
                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*uxi - 2.0d0*rr5*qix
                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*uyi - 2.0d0*rr5*qiy
                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*uzi - 2.0d0*rr5*qiz
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)*dscale(kk)
                     field(j,k) = field(j,k) + fkd(j)*dscale(kk)
                     fieldp(j,i) = fieldp(j,i) + fid(j)*pscale(kk)
                     fieldp(j,k) = fieldp(j,k) + fkd(j)*pscale(kk)
                  end do
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     find Poisson-Boltzmann reaction field due to permanent multipoles
c
      call pbempole
c
c     perform dynamic allocation of some local arrays
c
      allocate (fields(3,npole))
      allocate (fieldps(3,npole))
      allocate (udir(3,npole))
      allocate (udirp(3,npole))
      allocate (udirs(3,npole))
      allocate (udirps(3,npole))
      allocate (uold(3,npole))
      allocate (uoldp(3,npole))
      allocate (uolds(3,npole))
      allocate (uoldps(3,npole))
      allocate (indpole(3,n))
      allocate (inppole(3,n))
c
c     set vacuum induced dipoles to polarizability times direct field;
c     SCRF induced dipoles are polarizability times direct field
c     plus the reaction field due to permanent multipoles
c
      do i = 1, npole
         ii = ipole(i)
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            udirs(j,i) = polarity(i) * (field(j,i)+pbep(j,ii))
            udirps(j,i) = polarity(i) * (fieldp(j,i)+pbep(j,ii))
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
            uinds(j,i) = udirs(j,i)
            uinps(j,i) = udirps(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         eps = 100.0d0
c
c     get mutual induced dipole moments by an iterative method
c
         do while (.not. done)
            do i = 1, n
               do j = 1, 3
                  indpole(j,i) = 0.0d0
                  inppole(j,i) = 0.0d0
                  pbeuind(j,i) = 0.0d0
                  pbeuinp(j,i) = 0.0d0
               end do
            end do
            do i = 1, npole
               ii = ipole(i)
               do j = 1, 3
                  indpole(j,ii) = uinds(j,i)
                  inppole(j,ii) = uinps(j,i)
               end do
            end do
            call apbsinduce (indpole,pbeuind)
            call apbsnlinduce (inppole,pbeuinp)
c
c     compute intra-solute mutual polarization
c
            do i = 1, npole
               do j = 1, 3
                  field(j,i) = 0.0d0
                  fieldp(j,i) = 0.0d0
                  fields(j,i) = 0.0d0
                  fieldps(j,i) = 0.0d0
               end do
            end do
            do i = 1, npole
               ii = ipole(i)
               pdi = pdamp(i)
               pti = thole(i)
               duix = uind(1,i)
               duiy = uind(2,i)
               duiz = uind(3,i)
               puix = uinp(1,i)
               puiy = uinp(2,i)
               puiz = uinp(3,i)
               duixs = uinds(1,i)
               duiys = uinds(2,i)
               duizs = uinds(3,i)
               puixs = uinps(1,i)
               puiys = uinps(2,i)
               puizs = uinps(3,i)
               do j = 1, np11(ii)
                  dscale(ip11(j,ii)) = u1scale
               end do
               do j = 1, np12(ii)
                  dscale(ip12(j,ii)) = u2scale
               end do
               do j = 1, np13(ii)
                  dscale(ip13(j,ii)) = u3scale
               end do
               do j = 1, np14(ii)
                  dscale(ip14(j,ii)) = u4scale
               end do
               do k = i+1, npole
                  kk = ipole(k)
                  proceed = .true.
                  if (use_intra)
     &               call groups (proceed,fgrp,ii,kk,0,0,0,0)
                  if (proceed) then
                     xr = x(kk) - x(ii)
                     yr = y(kk) - y(ii)
                     zr = z(kk) - z(ii)
                     xr2 = xr * xr
                     yr2 = yr * yr
                     zr2 = zr * zr
                     r2 = xr2 + yr2 + zr2
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        dukx = uind(1,k)
                        duky = uind(2,k)
                        dukz = uind(3,k)
                        pukx = uinp(1,k)
                        puky = uinp(2,k)
                        pukz = uinp(3,k)
                        dukxs = uinds(1,k)
                        dukys = uinds(2,k)
                        dukzs = uinds(3,k)
                        pukxs = uinps(1,k)
                        pukys = uinps(2,k)
                        pukzs = uinps(3,k)
                        scale3 = dscale(kk)
                        scale5 = dscale(kk)
                        damp = pdi * pdamp(k)
                        if (damp .ne. 0.0d0) then
                           pgamma = min(pti,thole(k))
                           damp = -pgamma * (r/damp)**3
                           if (damp .gt. -50.0d0) then
                              expdamp = exp(damp)
                              scale3 = scale3 * (1.0d0-expdamp)
                              scale5 = scale5 * (1.0d0-(1.0d0-damp)
     &                                              *expdamp)
                           end if
                        end if
                        rr3 = scale3 / (r*r2)
                        rr5 = 3.0d0 * scale5 / (r*r2*r2)
                        duir = xr*duix + yr*duiy + zr*duiz
                        dukr = xr*dukx + yr*duky + zr*dukz
                        puir = xr*puix + yr*puiy + zr*puiz
                        pukr = xr*pukx + yr*puky + zr*pukz
                        duirs = xr*duixs + yr*duiys + zr*duizs
                        dukrs = xr*dukxs + yr*dukys + zr*dukzs
                        puirs = xr*puixs + yr*puiys + zr*puizs
                        pukrs = xr*pukxs + yr*pukys + zr*pukzs
                        fid(1) = -rr3*dukx + rr5*dukr*xr
                        fid(2) = -rr3*duky + rr5*dukr*yr
                        fid(3) = -rr3*dukz + rr5*dukr*zr
                        fkd(1) = -rr3*duix + rr5*duir*xr
                        fkd(2) = -rr3*duiy + rr5*duir*yr
                        fkd(3) = -rr3*duiz + rr5*duir*zr
                        fip(1) = -rr3*pukx + rr5*pukr*xr
                        fip(2) = -rr3*puky + rr5*pukr*yr
                        fip(3) = -rr3*pukz + rr5*pukr*zr
                        fkp(1) = -rr3*puix + rr5*puir*xr
                        fkp(2) = -rr3*puiy + rr5*puir*yr
                        fkp(3) = -rr3*puiz + rr5*puir*zr
                        fids(1) = -rr3*dukxs + rr5*dukrs*xr
                        fids(2) = -rr3*dukys + rr5*dukrs*yr
                        fids(3) = -rr3*dukzs + rr5*dukrs*zr
                        fkds(1) = -rr3*duixs + rr5*duirs*xr
                        fkds(2) = -rr3*duiys + rr5*duirs*yr
                        fkds(3) = -rr3*duizs + rr5*duirs*zr
                        fips(1) = -rr3*pukxs + rr5*pukrs*xr
                        fips(2) = -rr3*pukys + rr5*pukrs*yr
                        fips(3) = -rr3*pukzs + rr5*pukrs*zr
                        fkps(1) = -rr3*puixs + rr5*puirs*xr
                        fkps(2) = -rr3*puiys + rr5*puirs*yr
                        fkps(3) = -rr3*puizs + rr5*puirs*zr
                        do j = 1, 3
                           field(j,i) = field(j,i) + fid(j)
                           field(j,k) = field(j,k) + fkd(j)
                           fieldp(j,i) = fieldp(j,i) + fip(j)
                           fieldp(j,k) = fieldp(j,k) + fkp(j)
                           fields(j,i) = fields(j,i) + fids(j)
                           fields(j,k) = fields(j,k) + fkds(j)
                           fieldps(j,i) = fieldps(j,i) + fips(j)
                           fieldps(j,k) = fieldps(j,k) + fkps(j)
                        end do
                     end if
                  end if
               end do
c
c     reset interaction scaling coefficients for connected atoms
c
               do j = 1, np11(ii)
                  dscale(ip11(j,ii)) = 1.0d0
               end do
               do j = 1, np12(ii)
                  dscale(ip12(j,ii)) = 1.0d0
               end do
               do j = 1, np13(ii)
                  dscale(ip13(j,ii)) = 1.0d0
               end do
               do j = 1, np14(ii)
                  dscale(ip14(j,ii)) = 1.0d0
               end do
            end do
c
c     check to see if the mutual induced dipoles have converged
c
            iter = iter + 1
            epsold = eps
            epsd = 0.0d0
            epsp = 0.0d0
            epsds = 0.0d0
            epsps = 0.0d0
            do i = 1, npole
               ii = ipole(i)
               do j = 1, 3
                  uold(j,i) = uind(j,i)
                  uoldp(j,i) = uinp(j,i)
                  uolds(j,i) = uinds(j,i)
                  uoldps(j,i) = uinps(j,i)
                  uind(j,i) = udir(j,i) + polarity(i)*field(j,i)
                  uinp(j,i) = udirp(j,i) + polarity(i)*fieldp(j,i)
                  uinds(j,i) = udirs(j,i) + polarity(i)
     &                              *(fields(j,i)+pbeuind(j,ii))
                  uinps(j,i) = udirps(j,i) + polarity(i)
     &                              *(fieldps(j,i)+pbeuinp(j,ii))
                  uind(j,i) = uold(j,i) + polsor*(uind(j,i)-uold(j,i))
                  uinp(j,i) = uoldp(j,i) + polsor*(uinp(j,i)-uoldp(j,i))
                  uinds(j,i) = uolds(j,i) + polsor
     &                              *(uinds(j,i)-uolds(j,i))
                  uinps(j,i) = uoldps(j,i) + polsor
     &                              *(uinps(j,i)-uoldps(j,i))
                  epsd = epsd + (uind(j,i)-uold(j,i))**2
                  epsp = epsp + (uinp(j,i)-uoldp(j,i))**2
                  epsds = epsds + (uinds(j,i)-uolds(j,i))**2
                  epsps = epsps + (uinps(j,i)-uoldps(j,i))**2
               end do
            end do
            eps = max(epsd,epsp,epsds,epsps)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. maxiter)  done = .true.
         end do
c
c     terminate the calculation if dipoles failed to converge
c
         if (eps .gt. poleps) then
            write (iout,30)
   30       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      deallocate (field)
      deallocate (fieldp)
      deallocate (fields)
      deallocate (fieldps)
      deallocate (udir)
      deallocate (udirp)
      deallocate (udirs)
      deallocate (udirps)
      deallocate (uold)
      deallocate (uoldp)
      deallocate (uolds)
      deallocate (uoldps)
      deallocate (indpole)
      deallocate (inppole)
      return
      end
