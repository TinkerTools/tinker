c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  ewreg.i  --  exponential factors for regular Ewald sum  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     maxvec    maximum number of k-vectors per reciprocal axis
c
c     ejc       exponental factors for cosine along the j-axis
c     ejs       exponental factors for sine along the j-axis
c     ekc       exponental factors for cosine along the k-axis
c     eks       exponental factors for sine along the k-axis
c     elc       exponental factors for cosine along the l-axis
c     els       exponental factors for sine along the l-axis
c
c
      integer maxvec
      parameter (maxvec=15)
      real*8 ejc,ejs,ekc,eks,elc,els
      common /ewreg/ ejc(maxatm,0:maxvec),ejs(maxatm,0:maxvec),
     &               ekc(maxatm,-maxvec:maxvec),
     &               eks(maxatm,-maxvec:maxvec),
     &               elc(maxatm,-maxvec:maxvec),
     &               els(maxatm,-maxvec:maxvec)
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kewald  --  Ewald sum parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kewald" assigns both regular Ewald summation and particle
c     mesh Ewald parameters for a periodic box
c
c
      subroutine kewald
      implicit none
      include 'sizes.i'
      include 'boxes.i'
      include 'cutoff.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'pme.i'
      include 'potent.i'
      integer maxpower
      parameter (maxpower=50)
      integer i,k
      integer next,ntable
      integer ifft1,ifft2,ifft3
      integer jmax,kmax,lmax
      integer multi(maxpower)
      real*8 delta,eps
      logical use_pme
      character*20 keyword
      character*120 record
      character*120 string
      data multi  /   2,   3,   4,   5,   6,   8,   9,  10,  12,  15,
     &               16,  18,  20,  24,  25,  27,  30,  32,  36,  40,
     &               45,  48,  50,  54,  60,  64,  72,  75,  80,  81,
     &               90,  96, 100, 108, 120, 125, 128, 135, 144, 150,
     &              160, 162, 180, 192, 200, 216, 225, 240, 243, 250 /
c
c
c     decide between regular Ewald and particle mesh Ewald
c
      use_pme = .false.
      if (use_charge)  use_pme = .true.
c
c     set defaults for reciprocal fraction and boundary treatment
c
      frecip = 0.5d0
      tinfoil = .true.
c
c     estimate an optimal value for the Ewald coefficient
c
      eps = 1.0d-8
      call ewaldcof (aewald,ewaldcut,eps)
c
c     set defaults for PME B-spline order and charge grid size;
c     grid is a product of powers of 2, 3 and/or 5 for efficiency
c
      if (use_pme) then
         bsorder = 8
         delta = 1.0d-8
         nfft1 = 0
         nfft2 = 0
         nfft3 = 0
         ifft1 = int(1.5d0*xbox-delta) + 1
         ifft2 = int(1.5d0*ybox-delta) + 1
         ifft3 = int(1.5d0*zbox-delta) + 1
         do i = maxpower, 1, -1
            k = multi(i)
            if (k .le. maxfft) then
               if (k.ge.ifft1 .or. nfft1.eq.0)  nfft1 = k
               if (k.ge.ifft2 .or. nfft2.eq.0)  nfft2 = k
               if (k.ge.ifft3 .or. nfft3.eq.0)  nfft3 = k
            end if
         end do
      end if
c
c     search keywords for Ewald summation commands
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:12) .eq. 'EWALD-ALPHA ') then
            read (string,*,err=20,end=20)  aewald
         else if (keyword(1:15) .eq. 'EWALD-FRACTION ') then
            read (string,*,err=20,end=20)  frecip
         else if (keyword(1:15) .eq. 'EWALD-BOUNDARY ') then
            tinfoil = .false.
         else if (keyword(1:9) .eq. 'PME-GRID ') then
            nfft1 = 0
            nfft2 = 0
            nfft3 = 0
            read (string,*,err=10,end=10)  nfft1,nfft2,nfft3
   10       continue
            if (nfft2 .eq. 0)  nfft2 = nfft1
            if (nfft3 .eq. 0)  nfft3 = nfft1
         else if (keyword(1:10) .eq. 'PME-ORDER ') then
            read (string,*,err=20,end=20)  bsorder
         end if
   20    continue
      end do
c
c     check the B-spline order and charge grid dimension for PME
c
      if (use_pme) then
         if (bsorder .gt. maxorder) then
            write (iout,30)
   30       format (/,' KEWALD  --  B-Spline Order Too Large;',
     &                 ' Increase MAXORDER')
            call fatal
         end if
         if (max(nfft1,nfft2,nfft3) .gt. maxfft) then
            write (iout,40)
   40       format (/,' KEWALD  --  FFT Charge Grid Too Large;',
     &                 ' Increase MAXFFT')
            call fatal
         else if (nfft1.lt.ifft1 .or. nfft2.lt.ifft2
     &                .or. nfft3.lt.ifft3) then
            write (iout,50)
   50       format (/,' KEWALD  --  Warning, Small Charge Grid',
     &                 ' may give Poor Accuracy')
         end if
c
c     check the number of k-vectors for regular Ewald summation
c
      else
         jmax = int(frecip/recip(1,1))
         kmax = int(frecip/recip(2,2))
         lmax = int(frecip/recip(3,3))
         if (max(jmax,kmax,lmax) .gt. maxvec) then
            jmax = min(maxvec,jmax)
            kmax = min(maxvec,kmax)
            lmax = min(maxvec,lmax)
            write (iout,60)
   60       format (/,' KEWALD  --  Too many Reciprocal Space',
     &                 ' K-Vectors; Increase MAXVEC')
         end if
      end if
c
c     initialize the PME arrays that can be precomputed
c
      if (use_pme) then
         ntable = maxtable
         call moduli
         call fftsetup (nfft1,nfft2,nfft3,ntable,table)
      end if
c
c     print a message listing some of the Ewald parameters
c
      if (verbose) then
         if (use_pme) then
            write (iout,70)  aewald,nfft1,nfft2,nfft3,bsorder
   70       format (/,' Smooth Particle Mesh Ewald Parameters :',
     &              //,4x,'Ewald Coefficient',6x,'Charge Grid',
     &                 ' Dimensions',6x,'B-Spline Order',
     &              //,8x,f8.4,11x,3i6,12x,i6)
         else
            write (iout,80)  aewald,jmax,kmax,lmax,frecip
   80       format (/,' Regular Ewald Summation Parameters :',
     &              //,4x,'Ewald Coefficient',6x,'K-Vector',
     &                 ' Dimensions',5x,'Reciprocal Fraction',
     &              //,8x,f8.4,11x,3i5,14x,f8.4)
         end if
      end if
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
c     permanent atomic multipole moments to the electrostatic field
c     for use in finding the direct induced dipole moments via a
c     regular Ewald summation
c
c
      subroutine udirect1 (field)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'math.i'
      include 'mpole.i'
      include 'units.i'
      integer i,j,k,l,ii
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 expterm,cut
      real*8 term,fterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 qf,t1,t2
      real*8 ck,dk,qk
      real*8 q1,q2,q3
      real*8 ckr(maxatm)
      real*8 skr(maxatm)
      real*8 cjk(maxatm)
      real*8 sjk(maxatm)
      real*8 cm(maxatm)
      real*8 dm(3,maxatm)
      real*8 qm(9,maxatm)
      real*8 field(3,maxatm)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      term = -0.25d0 / aewald**2
      fterm = 8.0d0 * pi / volbox
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
      jmax = min(maxvec,int(frecip/recip(1,1)))
      kmax = min(maxvec,int(frecip/recip(2,2)))
      lmax = min(maxvec,int(frecip/recip(3,3)))
c
c     copy the multipole moments into local storage areas
c
      do i = 1, npole
         cm(i) = rpole(1,i)
         dm(1,i) = rpole(2,i)
         dm(2,i) = rpole(3,i)
         dm(3,i) = rpole(4,i)
         qm(1,i) = rpole(5,i)
         qm(2,i) = rpole(6,i)
         qm(3,i) = rpole(7,i)
         qm(4,i) = rpole(8,i)
         qm(5,i) = rpole(9,i)
         qm(6,i) = rpole(10,i)
         qm(7,i) = rpole(11,i)
         qm(8,i) = rpole(12,i)
         qm(9,i) = rpole(13,i)
      end do
c
c     calculate and store the exponential factors
c
      do i = 1, npole
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(i,0) = 1.0d0
         ejs(i,0) = 0.0d0
         ekc(i,0) = 1.0d0
         eks(i,0) = 0.0d0
         elc(i,0) = 1.0d0
         els(i,0) = 0.0d0
         ejc(i,1) = cos(xfr)
         ejs(i,1) = sin(xfr)
         ekc(i,1) = cos(yfr)
         eks(i,1) = sin(yfr)
         elc(i,1) = cos(zfr)
         els(i,1) = sin(zfr)
         ekc(i,-1) = ekc(i,1)
         eks(i,-1) = -eks(i,1)
         elc(i,-1) = elc(i,1)
         els(i,-1) = -els(i,1)
         do j = 2, jmax
            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
         end do
         do j = 2, kmax
            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
            ekc(i,-j) = ekc(i,j)
            eks(i,-j) = -eks(i,j)
         end do
         do j = 2, lmax
            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
            elc(i,-j) = elc(i,j)
            els(i,-j) = -els(i,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do i = 1, npole
               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  do i = 1, npole
                     ckr(i) = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr(i) = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                     ck = cm(i)
                     dk = h1*dm(1,i) + h2*dm(2,i) + h3*dm(3,i)
                     q1 = h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i)
                     q2 = h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i)
                     q3 = h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i)
                     qk = h1*q1 + h2*q2 + h3*q3
                     t1 = t1 + (ck-qk)*skr(i) + dk*ckr(i)
                     t2 = t2 + (ck-qk)*ckr(i) - dk*skr(i)
                  end do
                  expterm = fterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  do i = 1, npole
                     qf = expterm * (skr(i)*t2-ckr(i)*t1)
                     field(1,i) = field(1,i) + h1*qf
                     field(2,i) = field(2,i) + h2*qf
                     field(3,i) = field(3,i) + h3*qf
                  end do
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
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
c     induced atomic dipole moments to the electrostatic field for
c     use in iterative calculation of induced dipole moments via a
c     regular Ewald summation
c
c
      subroutine umutual1 (field,fieldp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'units.i'
      integer i,j,k,l,ii
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 expterm,cut
      real*8 term,fterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 duk,puk
      real*8 dqf,pqf
      real*8 t1,t2,t3,t4
      real*8 ckr(maxatm)
      real*8 skr(maxatm)
      real*8 cjk(maxatm)
      real*8 sjk(maxatm)
      real*8 field(3,maxatm)
      real*8 fieldp(3,maxatm)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      term = -0.25d0 / aewald**2
      fterm = 8.0d0 * pi / volbox
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
      jmax = min(maxvec,int(frecip/recip(1,1)))
      kmax = min(maxvec,int(frecip/recip(2,2)))
      lmax = min(maxvec,int(frecip/recip(3,3)))
c
c     calculate and store the exponential factors
c
      do i = 1, npole
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(i,0) = 1.0d0
         ejs(i,0) = 0.0d0
         ekc(i,0) = 1.0d0
         eks(i,0) = 0.0d0
         elc(i,0) = 1.0d0
         els(i,0) = 0.0d0
         ejc(i,1) = cos(xfr)
         ejs(i,1) = sin(xfr)
         ekc(i,1) = cos(yfr)
         eks(i,1) = sin(yfr)
         elc(i,1) = cos(zfr)
         els(i,1) = sin(zfr)
         ekc(i,-1) = ekc(i,1)
         eks(i,-1) = -eks(i,1)
         elc(i,-1) = elc(i,1)
         els(i,-1) = -els(i,1)
         do j = 2, jmax
            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
         end do
         do j = 2, kmax
            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
            ekc(i,-j) = ekc(i,j)
            eks(i,-j) = -eks(i,j)
         end do
         do j = 2, lmax
            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
            elc(i,-j) = elc(i,j)
            els(i,-j) = -els(i,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do i = 1, npole
               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  t3 = 0.0d0
                  t4 = 0.0d0
                  do i = 1, npole
                     ckr(i) = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr(i) = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                     duk = h1*uind(1,i) + h2*uind(2,i) + h3*uind(3,i)
                     puk = h1*uinp(1,i) + h2*uinp(2,i) + h3*uinp(3,i)
                     t1 = t1 + duk*ckr(i)
                     t2 = t2 - duk*skr(i)
                     t3 = t3 + puk*ckr(i)
                     t4 = t4 - puk*skr(i)
                  end do
                  expterm = fterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  do i = 1, npole
                     dqf = expterm * (skr(i)*t2-ckr(i)*t1)
                     pqf = expterm * (skr(i)*t4-ckr(i)*t3)
                     field(1,i) = field(1,i) + h1*dqf
                     field(2,i) = field(2,i) + h2*dqf
                     field(3,i) = field(3,i) + h3*dqf
                     fieldp(1,i) = fieldp(1,i) + h1*pqf
                     fieldp(2,i) = fieldp(2,i) + h2*pqf
                     fieldp(3,i) = fieldp(3,i) + h3*pqf
                  end do
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine erecip  --  multipole Ewald reciprocal energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "erecip" evaluates the reciprocal space portion of the regular
c     Ewald summation energy due to atomic multipole interactions
c     and dipole polarizability
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine erecip
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'units.i'
      integer i,j,k,l,ii
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 e,ei,f,term,cut
      real*8 eterm,expterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 t1,t2,t3,t4
      real*8 ck,dk,qk,uk
      real*8 q1,q2,q3
      real*8 ckr,skr
      real*8 cm(maxatm)
      real*8 dm(3,maxatm)
      real*8 qm(9,maxatm)
      real*8 um(3,maxatm)
      real*8 cjk(maxatm)
      real*8 sjk(maxatm)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
      term = -0.25d0 / aewald**2
      eterm = 4.0d0 * pi * f / volbox
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
      jmax = min(maxvec,int(frecip/recip(1,1)))
      kmax = min(maxvec,int(frecip/recip(2,2)))
      lmax = min(maxvec,int(frecip/recip(3,3)))
c
c     copy the multipole moments into local storage areas
c
      do i = 1, npole
         cm(i) = rpole(1,i)
         dm(1,i) = rpole(2,i)
         dm(2,i) = rpole(3,i)
         dm(3,i) = rpole(4,i)
         qm(1,i) = rpole(5,i)
         qm(2,i) = rpole(6,i)
         qm(3,i) = rpole(7,i)
         qm(4,i) = rpole(8,i)
         qm(5,i) = rpole(9,i)
         qm(6,i) = rpole(10,i)
         qm(7,i) = rpole(11,i)
         qm(8,i) = rpole(12,i)
         qm(9,i) = rpole(13,i)
         um(1,i) = uind(1,i)
         um(2,i) = uind(2,i)
         um(3,i) = uind(3,i)
      end do
c
c     calculate and store the exponential factors
c
      do i = 1, npole
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(i,0) = 1.0d0
         ejs(i,0) = 0.0d0
         ekc(i,0) = 1.0d0
         eks(i,0) = 0.0d0
         elc(i,0) = 1.0d0
         els(i,0) = 0.0d0
         ejc(i,1) = cos(xfr)
         ejs(i,1) = sin(xfr)
         ekc(i,1) = cos(yfr)
         eks(i,1) = sin(yfr)
         elc(i,1) = cos(zfr)
         els(i,1) = sin(zfr)
         ekc(i,-1) = ekc(i,1)
         eks(i,-1) = -eks(i,1)
         elc(i,-1) = elc(i,1)
         els(i,-1) = -els(i,1)
         do j = 2, jmax
            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
         end do
         do j = 2, kmax
            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
            ekc(i,-j) = ekc(i,j)
            eks(i,-j) = -eks(i,j)
         end do
         do j = 2, lmax
            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
            elc(i,-j) = elc(i,j)
            els(i,-j) = -els(i,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do i = 1, npole
               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  t3 = 0.0d0
                  t4 = 0.0d0
                  do i = 1, npole
                     ckr = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                     ck = cm(i)
                     dk = h1*dm(1,i) + h2*dm(2,i) + h3*dm(3,i)
                     q1 = h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i)
                     q2 = h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i)
                     q3 = h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i)
                     qk = h1*q1 + h2*q2 + h3*q3
                     uk = h1*um(1,i) + h2*um(2,i) + h3*um(3,i)
                     t1 = t1 + (ck-qk)*skr + dk*ckr
                     t2 = t2 + (ck-qk)*ckr - dk*skr
                     t3 = t3 + uk*ckr
                     t4 = t4 - uk*skr
                  end do
                  expterm = eterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  e = expterm * (t1*t1+t2*t2)
                  ei = expterm * (t1*t3+t2*t4)
                  em = em + e
                  ep = ep + ei
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine erecip1  --  mpole Ewald recip energy & derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erecip1" evaluates the reciprocal space portion of the regular
c     Ewald summation energy and gradient due to atomic multipole
c     interactions and dipole polarizability
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine erecip1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'units.i'
      include 'virial.i'
      integer i,j,k,l
      integer ii,m1,m2
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 e,ei,etot,f,cut
      real*8 expterm,term,eterm
      real*8 uterm,vterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 ck,dk,qk,uk
      real*8 q1,q2,q3
      real*8 t1,t2,t3,t4
      real*8 ukp,t3p,t4p
      real*8 de,det1,det2
      real*8 dei,det1i,det2i
      real*8 wterm(3,3)
      real*8 t5(3,3),t6(3,3)
      real*8 t5u(3,3),t5p(3,3)
      real*8 t6u(3,3),t6p(3,3)
      real*8 qt(3,3),dt(3,3)
      real*8 dtu(3,3),dtp(3,3)
      real*8 ckr(maxatm),skr(maxatm)
      real*8 cjk(maxatm),sjk(maxatm)
      real*8 s1(maxatm),s2(maxatm)
      real*8 s3(maxatm),s4(maxatm)
      real*8 s3p(maxatm),s4p(maxatm)
      real*8 cm(maxatm),dm(3,maxatm)
      real*8 qm(9,maxatm),um(3,maxatm)
      real*8 trq(3,maxatm),trqi(3,maxatm)
      real*8 dkx(maxatm),qkx(maxatm)
      real*8 dky(maxatm),qky(maxatm)
      real*8 dkz(maxatm),qkz(maxatm)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
      term = -0.25d0 / aewald**2
      eterm = 4.0d0 * pi * f / volbox
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
      jmax = min(maxvec,int(frecip/recip(1,1)))
      kmax = min(maxvec,int(frecip/recip(2,2)))
      lmax = min(maxvec,int(frecip/recip(3,3)))
c
c     copy the multipole moments into local storage areas
c
      do i = 1, npole
         cm(i) = rpole(1,i)
         dm(1,i) = rpole(2,i)
         dm(2,i) = rpole(3,i)
         dm(3,i) = rpole(4,i)
         qm(1,i) = rpole(5,i)
         qm(2,i) = rpole(6,i)
         qm(3,i) = rpole(7,i)
         qm(4,i) = rpole(8,i)
         qm(5,i) = rpole(9,i)
         qm(6,i) = rpole(10,i)
         qm(7,i) = rpole(11,i)
         qm(8,i) = rpole(12,i)
         qm(9,i) = rpole(13,i)
         um(1,i) = uind(1,i)
         um(2,i) = uind(2,i)
         um(3,i) = uind(3,i)
      end do
c
c     zero out local accumulation arrays for derivatives
c
      do i = 1, n
         trq(1,i) = 0.0d0
         trq(2,i) = 0.0d0
         trq(3,i) = 0.0d0
         trqi(1,i) = 0.0d0
         trqi(2,i) = 0.0d0
         trqi(3,i) = 0.0d0
      end do
c
c     calculate and store the exponential factors
c
      do i = 1, npole
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(i,0) = 1.0d0
         ejs(i,0) = 0.0d0
         ekc(i,0) = 1.0d0
         eks(i,0) = 0.0d0
         elc(i,0) = 1.0d0
         els(i,0) = 0.0d0
         ejc(i,1) = cos(xfr)
         ejs(i,1) = sin(xfr)
         ekc(i,1) = cos(yfr)
         eks(i,1) = sin(yfr)
         elc(i,1) = cos(zfr)
         els(i,1) = sin(zfr)
         ekc(i,-1) = ekc(i,1)
         eks(i,-1) = -eks(i,1)
         elc(i,-1) = elc(i,1)
         els(i,-1) = -els(i,1)
         do j = 2, jmax
            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
         end do
         do j = 2, kmax
            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
            ekc(i,-j) = ekc(i,j)
            eks(i,-j) = -eks(i,j)
         end do
         do j = 2, lmax
            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
            elc(i,-j) = elc(i,j)
            els(i,-j) = -els(i,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do i = 1, npole
               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  t3 = 0.0d0
                  t4 = 0.0d0
                  t3p = 0.0d0
                  t4p = 0.0d0
                  do m2 = 1, 3
                     do m1 = 1, 3
                        t5(m1,m2) = 0.0d0
                        t5u(m1,m2) = 0.0d0
                        t5p(m1,m2) = 0.0d0
                        t6(m1,m2) = 0.0d0
                        t6u(m1,m2) = 0.0d0
                        t6p(m1,m2) = 0.0d0
                     end do
                  end do
                  do i = 1, npole
                     ckr(i) = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr(i) = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                     ck = cm(i)
                     dk = h1*dm(1,i) + h2*dm(2,i) + h3*dm(3,i)
                     dkx(i) = h3*dm(2,i) - h2*dm(3,i)
                     dky(i) = h1*dm(3,i) - h3*dm(1,i)
                     dkz(i) = h2*dm(1,i) - h1*dm(2,i)
                     q1 = h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i)
                     q2 = h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i)
                     q3 = h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i)
                     qk = h1*q1 + h2*q2 + h3*q3
                     qkx(i) = h3*q2 - h2*q3
                     qky(i) = h1*q3 - h3*q1
                     qkz(i) = h2*q1 - h1*q2
                     uk = h1*uind(1,i) + h2*uind(2,i) + h3*uind(3,i)
                     ukp = h1*uinp(1,i) + h2*uinp(2,i) + h3*uinp(3,i)
                     s1(i) = (ck-qk)*skr(i) + dk*ckr(i)
                     s2(i) = (ck-qk)*ckr(i) - dk*skr(i)
                     s3(i) = uk * ckr(i)
                     s4(i) = -uk * skr(i)
                     s3p(i) = ukp * ckr(i)
                     s4p(i) = -ukp * skr(i)
                     t1 = t1 + s1(i)
                     t2 = t2 + s2(i)
                     t3 = t3 + s3(i)
                     t4 = t4 + s4(i)
                     t3p = t3p + s3p(i)
                     t4p = t4p + s4p(i)
c
c     terms needed for subsequent virial tensor calculation
c
                     qt(1,1) = h1*(h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i))
                     qt(2,1) = h1*(h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i))
                     qt(3,1) = h1*(h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i))
                     qt(1,2) = h2*(h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i))
                     qt(2,2) = h2*(h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i))
                     qt(3,2) = h2*(h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i))
                     qt(1,3) = h3*(h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i))
                     qt(2,3) = h3*(h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i))
                     qt(3,3) = h3*(h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i))
                     dt(1,1) = h1 * dm(1,i)
                     dt(2,1) = h1 * dm(2,i)
                     dt(3,1) = h1 * dm(3,i)
                     dt(1,2) = h2 * dm(1,i)
                     dt(2,2) = h2 * dm(2,i)
                     dt(3,2) = h2 * dm(3,i)
                     dt(1,3) = h3 * dm(1,i)
                     dt(2,3) = h3 * dm(2,i)
                     dt(3,3) = h3 * dm(3,i)
                     dtu(1,1) = h1 * uind(1,i)
                     dtu(2,1) = h1 * uind(2,i)
                     dtu(3,1) = h1 * uind(3,i)
                     dtu(1,2) = h2 * uind(1,i)
                     dtu(2,2) = h2 * uind(2,i)
                     dtu(3,2) = h2 * uind(3,i)
                     dtu(1,3) = h3 * uind(1,i)
                     dtu(2,3) = h3 * uind(2,i)
                     dtu(3,3) = h3 * uind(3,i)
                     dtp(1,1) = h1 * uinp(1,i)
                     dtp(2,1) = h1 * uinp(2,i)
                     dtp(3,1) = h1 * uinp(3,i)
                     dtp(1,2) = h2 * uinp(1,i)
                     dtp(2,2) = h2 * uinp(2,i)
                     dtp(3,2) = h2 * uinp(3,i)
                     dtp(1,3) = h3 * uinp(1,i)
                     dtp(2,3) = h3 * uinp(2,i)
                     dtp(3,3) = h3 * uinp(3,i)
                     do m2 = 1, 3
                        do m1 = 1, 3
                           t5(m1,m2) = t5(m1,m2) - dt(m1,m2)*ckr(i)
     &                                    + 2.0d0*qt(m1,m2)*skr(i)
                           t5u(m1,m2) = t5u(m1,m2) - dtu(m1,m2)*ckr(i)
                           t5p(m1,m2) = t5p(m1,m2) - dtp(m1,m2)*ckr(i)
                           t6(m1,m2) = t6(m1,m2) + dt(m1,m2)*skr(i)
     &                                    + 2.0d0*qt(m1,m2)*ckr(i)
                           t6u(m1,m2) = t6u(m1,m2) + dtu(m1,m2)*skr(i)
                           t6p(m1,m2) = t6p(m1,m2) + dtp(m1,m2)*skr(i)
                        end do
                     end do
                  end do
c
c     get the energy contributions for current reciprocal vector
c
                  expterm = eterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  e = expterm * (t1*t1+t2*t2)
                  ei = expterm * (t1*t3+t2*t4)
                  etot = e + ei
                  em = em + e
                  ep = ep + ei
c
c     get the virial contributions for current reciprocal vector
c
                  uterm = expterm * (t1*(t1+t3+t3p) + t3*t3p
     &                                 + t2*(t2+t4+t4p) + t4*t4p)
                  do m2 = 1, 3
                     do m1 = 1, 3
                        wterm(m1,m2) = 2.0d0 * expterm
     &                     * (t1*t5(m1,m2) + t2*t6(m1,m2)
     &                        + 0.5d0*(t1*(t5u(m1,m2)+t5p(m1,m2))
     &                        + t2*(t6u(m1,m2)+t6p(m1,m2))
     &                        + (t3+t3p)*t5(m1,m2)
     &                        + t3*t5p(m1,m2) + t3p*t5u(m1,m2)
     &                        + (t4+t4p)*t6(m1,m2)
     &                        + t4*t6p(m1,m2) + t4p*t6u(m1,m2)))
                     end do
                  end do
                  if (poltyp .eq. 'DIRECT') then
                     uterm = uterm - expterm*(t3*t3p+t4*t4p)
                     do m2 = 1, 3
                        do m1 = 1, 3
                           wterm(m1,m2) = wterm(m1,m2)
     &                        - expterm*(t3*t5p(m1,m2)+t3p*t5u(m1,m2)
     &                                  +t4*t6p(m1,m2)+t4p*t6u(m1,m2))
                        end do
                     end do
                  end if
                  wterm(2,1) = 0.5d0 * (wterm(2,1)+wterm(1,2))
                  wterm(3,1) = 0.5d0 * (wterm(3,1)+wterm(1,3))
                  wterm(3,2) = 0.5d0 * (wterm(3,2)+wterm(2,3))
                  wterm(1,2) = wterm(2,1)
                  wterm(1,3) = wterm(3,1)
                  wterm(2,3) = wterm(3,2)
                  vterm = 2.0d0 * uterm * (1.0d0-term*hsq) / hsq
                  vir(1,1) = vir(1,1) + h1*h1*vterm + wterm(1,1) - uterm
                  vir(2,1) = vir(2,1) + h2*h1*vterm + wterm(2,1)
                  vir(3,1) = vir(3,1) + h3*h1*vterm + wterm(3,1)
                  vir(1,2) = vir(1,2) + h1*h2*vterm + wterm(1,2)
                  vir(2,2) = vir(2,2) + h2*h2*vterm + wterm(2,2) - uterm
                  vir(3,2) = vir(3,2) + h3*h2*vterm + wterm(3,2)
                  vir(1,3) = vir(1,3) + h1*h3*vterm + wterm(1,3)
                  vir(2,3) = vir(2,3) + h2*h3*vterm + wterm(2,3)
                  vir(3,3) = vir(3,3) + h3*h3*vterm + wterm(3,3) - uterm
c
c     get the force contributions for current reciprocal vector
c
                  expterm = 2.0d0 * expterm
                  do i = 1, npole
                     ii = ipole(i)
                     de = expterm * (s2(i)*t1-s1(i)*t2)
                     dei = 0.5d0 * expterm * ((s4(i)+s4p(i))*t1
     &                                       -(s3(i)+s3p(i))*t2
     &                                +s2(i)*(t3+t3p)-s1(i)*(t4+t4p))
                     if (poltyp .eq. 'MUTUAL') then
                         dei = dei + 0.5d0 * expterm
     &                            * (s4p(i)*t3+s4(i)*t3p
     &                              -s3p(i)*t4-s3(i)*t4p)
                     end if
                     det1 = expterm * (skr(i)*t2-ckr(i)*t1)
                     det2 = 2.0d0 * expterm * (ckr(i)*t2+skr(i)*t1)
                     det1i = 0.5d0 * expterm * (skr(i)*(t4+t4p)
     &                                         -ckr(i)*(t3+t3p))
                     det2i = expterm * (ckr(i)*(t4+t4p)+skr(i)*(t3+t3p))
                     dem(1,ii) = dem(1,ii) + h1*de
                     dem(2,ii) = dem(2,ii) + h2*de
                     dem(3,ii) = dem(3,ii) + h3*de
                     dep(1,ii) = dep(1,ii) + h1*dei
                     dep(2,ii) = dep(2,ii) + h2*dei
                     dep(3,ii) = dep(3,ii) + h3*dei
                     trq(1,ii) = trq(1,ii) + dkx(i)*det1 + qkx(i)*det2
                     trq(2,ii) = trq(2,ii) + dky(i)*det1 + qky(i)*det2
                     trq(3,ii) = trq(3,ii) + dkz(i)*det1 + qkz(i)*det2
                     trqi(1,ii) = trqi(1,ii) + dkx(i)*det1i
     &                               + qkx(i)*det2i
                     trqi(2,ii) = trqi(2,ii) + dky(i)*det1i
     &                               + qky(i)*det2i
                     trqi(3,ii) = trqi(3,ii) + dkz(i)*det1i
     &                               + qkz(i)*det2i
                  end do
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
c
c     convert the torques to forces and increment the totals
c
      call torque (trq,dem)
      call torque (trqi,dep)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine torque  --  convert multipole torque to force  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "torque" takes the torque values on sites defined by local
c     coordinate frames and distributes them to convert to forces
c     on the original sites and sites specifying the local frames
c
c     literature reference:
c
c     A. J. Stone and M. Alderton, "Distributed Multipole Analysis -
c     Methods and Applications", Molecular Physics, 56, 1047-1064 (1985)
c
c
      subroutine torque (trq,derivs)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i,j,ia,ib,ic
      real*8 usiz,vsiz,wsiz
      real*8 upsiz,vpsiz
      real*8 dotdu,dotdv
      real*8 dphidu,dphidv,dphidw
      real*8 c,s,uvdis,vudis,du,dv
      real*8 u(3),v(3),w(3)
      real*8 up(3),vp(3),diff(3)
      real*8 trq(3,maxatm)
      real*8 derivs(3,maxatm)
c
c
c     coordinate frame motion described by rotation about u, v and w
c
      do i = 1, npole
         ia = zaxis(i)
         ib = ipole(i)
         ic = xaxis(i)
         u(1) = x(ia) - x(ib)
         u(2) = y(ia) - y(ib)
         u(3) = z(ia) - z(ib)
         usiz = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))
         v(1) = x(ic) - x(ib)
         v(2) = y(ic) - y(ib)
         v(3) = z(ic) - z(ib)
         vsiz = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
         w(1) = u(2)*v(3) - u(3)*v(2)
         w(2) = u(3)*v(1) - u(1)*v(3)
         w(3) = u(1)*v(2) - u(2)*v(1)
         wsiz = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
         dotdu = 0.0d0
         dotdv = 0.0d0
         do j = 1, 3
            u(j) = u(j) / usiz
            v(j) = v(j) / vsiz
            w(j) = w(j) / wsiz
            diff(j) = v(j) - u(j)
            dotdu = dotdu + u(j)*diff(j)
            dotdv = dotdv + v(j)*diff(j)
         end do
c
c     get perpendiculars to u,v to get direction of motion of
c     u or v due to rotation about the cross product vector w
c
         upsiz = 0.0d0
         vpsiz = 0.0d0
         do j = 1, 3
            up(j) = diff(j) - dotdu*u(j)
            vp(j) = diff(j) - dotdv*v(j)
            upsiz = upsiz + up(j)*up(j)
            vpsiz = vpsiz + vp(j)*vp(j)
         end do
         upsiz = sqrt(upsiz)
         vpsiz = sqrt(vpsiz)
         do j = 1, 3
            up(j) = up(j) / upsiz
            vp(j) = vp(j) / vpsiz
         end do
c
c     negative of dot product of torque with unit vectors along u, v
c     and w give result of infinitesmal rotation along these vectors
c     i.e. dphi/dtheta = dot product, where dphi is work, dtheta is
c     angle; then dphi/dtheta is torque and the dot product is torque
c     component along unit vector
c
         dphidu = -trq(1,ib)*u(1) - trq(2,ib)*u(2) - trq(3,ib)*u(3)
         dphidv = -trq(1,ib)*v(1) - trq(2,ib)*v(2) - trq(3,ib)*v(3)
         dphidw = -trq(1,ib)*w(1) - trq(2,ib)*w(2) - trq(3,ib)*w(3)
c
c     get the projected distances between the vectors
c
         c = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
         s = sqrt(1.0d0 - c*c)
         uvdis = usiz * s
         vudis = vsiz * s
c
c     force distribution for the bisector local coordinate method
c
         if (polaxe(i) .eq. 'Bisector') then
            do j = 1, 3
               du = -w(j)*dphidv/uvdis + up(j)*dphidw/(2.0d0*usiz)
               dv = w(j)*dphidu/vudis + vp(j)*dphidw/(2.0d0*vsiz)
               derivs(j,ia) = derivs(j,ia) + du
               derivs(j,ic) = derivs(j,ic) + dv
               derivs(j,ib) = derivs(j,ib) - dv - du
            end do
c
c     force distribution for the Z-then-X local coordinate method
c
         else if (polaxe(i) .eq. 'Z-then-X') then
            do j = 1, 3
               du = -w(j)*dphidv/uvdis + up(j)*dphidw/usiz
               dv = w(j)*dphidu/vudis
               derivs(j,ia) = derivs(j,ia) + du
               derivs(j,ic) = derivs(j,ic) + dv
               derivs(j,ib) = derivs(j,ib) - dv - du
            end do
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine erecip3  --  mpole Ewald reciprocal analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "erecip3" evaluates the reciprocal space portion of the regular
c     Ewald summation energy due to atomic multipole interactions
c     and dipole polarizability, and prints information about the
c     energy over the reciprocal lattice vectors
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine erecip3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'units.i'
      integer i,j,k,l,ii
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 e,ei,f,term,cut
      real*8 eterm,expterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3
      real*8 hsq,hleng
      real*8 t1,t2,t3,t4
      real*8 ck,dk,qk,uk
      real*8 q1,q2,q3
      real*8 ckr,skr
      real*8 cm(maxatm)
      real*8 dm(3,maxatm)
      real*8 qm(9,maxatm)
      real*8 um(3,maxatm)
      real*8 cjk(maxatm)
      real*8 sjk(maxatm)
      logical header,huge
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      header = .true.
      f = electric / dielec
      term = -0.25d0 / aewald**2
      eterm = 4.0d0 * pi * f / volbox
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
      jmax = min(maxvec,int(frecip/recip(1,1)))
      kmax = min(maxvec,int(frecip/recip(2,2)))
      lmax = min(maxvec,int(frecip/recip(3,3)))
c
c     copy the multipole moments into local storage areas
c
      do i = 1, npole
         cm(i) = rpole(1,i)
         dm(1,i) = rpole(2,i)
         dm(2,i) = rpole(3,i)
         dm(3,i) = rpole(4,i)
         qm(1,i) = rpole(5,i)
         qm(2,i) = rpole(6,i)
         qm(3,i) = rpole(7,i)
         qm(4,i) = rpole(8,i)
         qm(5,i) = rpole(9,i)
         qm(6,i) = rpole(10,i)
         qm(7,i) = rpole(11,i)
         qm(8,i) = rpole(12,i)
         qm(9,i) = rpole(13,i)
         um(1,i) = uind(1,i)
         um(2,i) = uind(2,i)
         um(3,i) = uind(3,i)
      end do
c
c     calculate and store the exponential factors
c
      do i = 1, npole
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(i,0) = 1.0d0
         ejs(i,0) = 0.0d0
         ekc(i,0) = 1.0d0
         eks(i,0) = 0.0d0
         elc(i,0) = 1.0d0
         els(i,0) = 0.0d0
         ejc(i,1) = cos(xfr)
         ejs(i,1) = sin(xfr)
         ekc(i,1) = cos(yfr)
         eks(i,1) = sin(yfr)
         elc(i,1) = cos(zfr)
         els(i,1) = sin(zfr)
         ekc(i,-1) = ekc(i,1)
         eks(i,-1) = -eks(i,1)
         elc(i,-1) = elc(i,1)
         els(i,-1) = -els(i,1)
         do j = 2, jmax
            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
         end do
         do j = 2, kmax
            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
            ekc(i,-j) = ekc(i,j)
            eks(i,-j) = -eks(i,j)
         end do
         do j = 2, lmax
            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
            elc(i,-j) = elc(i,j)
            els(i,-j) = -els(i,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do i = 1, npole
               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  t3 = 0.0d0
                  t4 = 0.0d0
                  do i = 1, npole
                     ckr = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
c
c     alternative phase cosine and sine is slower but avoids storage
c
c                    ii = ipole(i)
c                    phi = x(i)*h1 + y(i)*h2 + z(i)*h3
c                    ckr = cos(phi)
c                    skr = sin(phi)
                     ck = cm(i)
                     dk = h1*dm(1,i) + h2*dm(2,i) + h3*dm(3,i)
                     q1 = h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i)
                     q2 = h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i)
                     q3 = h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i)
                     qk = h1*q1 + h2*q2 + h3*q3
                     uk = h1*um(1,i) + h2*um(2,i) + h3*um(3,i)
                     t1 = t1 + (ck-qk)*skr + dk*ckr
                     t2 = t2 + (ck-qk)*ckr - dk*skr
                     t3 = t3 + uk*ckr
                     t4 = t4 - uk*skr
                  end do
                  expterm = eterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  e = expterm * (t1*t1+t2*t2)
                  ei = expterm * (t1*t3+t2*t4)
                  nem = nem + 1
                  nep = nep + 1
                  em = em + e
                  ep = ep + ei
c
c     print a message if the energy of this interaction is large
c
                  huge = (max(abs(e),abs(ei)) .gt. 10.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Reciprocal Space Multipole',
     &                             ' and Polarization Terms :',
     &                          //,' Type',15x,'K-Vector',9x,'Fraction',
     &                             6x,'Energies (MPole, Polar)',/)
                     end if
                     hleng = hsq / (4.0d0*pi*pi)
                     write (iout,20)  j,k,l,hleng,e,ei
   20                format (' M-Pole',7x,3i5,7x,f8.4,4x,2f12.4)
                  end if
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
      return
      end
