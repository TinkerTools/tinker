c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung & Jay W. Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program testeost  --  test orthogonal space tempering    ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "testeost" checks helper routines used by the orthogonal
c     space tempering implementation against deterministic inputs
c     and expected outputs
c
c
      program testeost
      use bath
      use units
      implicit none
      integer nfail
c
c
c     initialize test state
c
      nfail = 0
      kelvin = 300.0d0
c
c     run individual tests
c
      call testindex (nfail)
      call testresize (nfail)
      call testbuildindex (nfail)
      call testensure (nfail)
      call testgkernels (nfail)
      call testfkernel (nfail)
      call testefkernel (nfail)
      call testmapsub (nfail)
      call testmeta (nfail)
c
c     print final result
c
      if (nfail .eq. 0) then
         write (*,10)
   10    format (/,' TESTEOST  --  All Tests Passed')
      else
         write (*,20)  nfail
   20    format (/,' TESTEOST  -- ',i6,' Tests Failed')
         stop 1
      end if
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine testindex  --  test packed index conversion   ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine testindex (nfail)
      implicit none
      integer nfail
      integer i,j,k
c
c
c     test ij_to_k and k_to_ij
c
      call ij_to_k (3,4,7,k)
      call checkint ('ij_to_k input i=3 j=4 nrow=7',
     &               k,24,nfail)
      call k_to_ij (24,7,i,j)
      call checkint ('k_to_ij input k=24 nrow=7 i',
     &               i,3,nfail)
      call checkint ('k_to_ij input k=24 nrow=7 j',
     &               j,4,nfail)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine testresize  --  test histogram resizing       ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine testresize (nfail)
      use ost
      implicit none
      integer i
      integer nfail
c
c
c     set small history arrays and resize them
c
      call resetost (5,5,2)
      nosthist = 2
      do i = 1, 2
         osthist(i) = 10 + i
         ostnext(i) = i - 1
         ostlhist(i) = 0.25d0 * dble(i)
         ostfhist(i) = -3.0d0 + 2.0d0*dble(i)
         osthhist(i) = 1.0d0 + dble(i)
         ostwlhist(i) = 0.25d0
         ostwfhist(i) = 1.0d0
      end do
      call resizeosthist
c
c     check old data are preserved and new slots are initialized
c
      call checkint ('resizeosthist size input size=2',
     &               sizeosthist,4,nfail)
      do i = 1, 4
         if (i .le. 2) then
            call checkint ('resizeosthist preserve osthist',
     &                     osthist(i),10+i,nfail)
            call checkint ('resizeosthist preserve ostnext',
     &                     ostnext(i),i-1,nfail)
            call checkreal ('resizeosthist preserve ostlhist',
     &                      ostlhist(i),0.25d0*dble(i),
     &                      1.0d-12,nfail)
            call checkreal ('resizeosthist preserve ostfhist',
     &                      ostfhist(i),-3.0d0+2.0d0*dble(i),
     &                      1.0d-12,nfail)
            call checkreal ('resizeosthist preserve osthhist',
     &                      osthhist(i),1.0d0+dble(i),
     &                      1.0d-12,nfail)
            call checkreal ('resizeosthist preserve ostwlhist',
     &                      ostwlhist(i),0.25d0,1.0d-12,nfail)
            call checkreal ('resizeosthist preserve ostwfhist',
     &                      ostwfhist(i),1.0d0,1.0d-12,nfail)
         else
            call checkint ('resizeosthist init osthist',
     &                     osthist(i),0,nfail)
            call checkint ('resizeosthist init ostnext',
     &                     ostnext(i),0,nfail)
            call checkreal ('resizeosthist init ostlhist',
     &                      ostlhist(i),0.0d0,1.0d-12,nfail)
            call checkreal ('resizeosthist init ostfhist',
     &                      ostfhist(i),0.0d0,1.0d-12,nfail)
            call checkreal ('resizeosthist init osthhist',
     &                      osthhist(i),0.0d0,1.0d-12,nfail)
            call checkreal ('resizeosthist init ostwlhist',
     &                      ostwlhist(i),0.0d0,1.0d-12,nfail)
            call checkreal ('resizeosthist init ostwfhist',
     &                      ostwfhist(i),0.0d0,1.0d-12,nfail)
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine testbuildindex -- test histogram index build  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine testbuildindex (nfail)
      use ost
      implicit none
      integer nfail
      integer k,ilmda,iflmda
c
c
c     put two gaussians in the same bin and one in another bin
c
      call resetost (5,5,3)
      nosthist = 3
      call sethist (1,0.50d0,0.0d0,1.0d0,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,2.0d0,wlmda,wflmda)
      call sethist (3,0.75d0,1.0d0,3.0d0,wlmda,wflmda)
      call buildostindex
c
c     expected linked list head is newest entry in the shared bin
c
      ilmda = 3
      iflmda = 3
      call ij_to_k (ilmda,iflmda,nlmda,k)
      call checkint ('buildostindex packed shared bin',
     &               osthist(1),k,nfail)
      call checkint ('buildostindex head shared bin',
     &               osthead(ilmda,iflmda),2,nfail)
      call checkint ('buildostindex next shared bin',
     &               ostnext(2),1,nfail)
      call checkint ('buildostindex tail shared bin',
     &               ostnext(1),0,nfail)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine testensure -- test flambda grid expansion     ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine testensure (nfail)
      use ost
      implicit none
      integer nfail
      integer k
      integer ilmda,iflmda
c
c
c     high-side expansion preserves fli0 and old gkernel values
c
      call resetost (3,5,1)
      gkernel(2,3) = 7.0d0
      call ensureflambda (2000.0d0)
      call checkint ('ensureflambda high nflmda',
     &               nflmda,3005,nfail)
      call checkint ('ensureflambda high fli0',
     &               fli0,3,nfail)
      call checkreal ('ensureflambda high preserve gkernel',
     &                gkernel(2,3),7.0d0,1.0d-12,nfail)
c
c     low-side expansion shifts fli0 and old gkernel values
c
      call resetost (3,5,1)
      gkernel(2,3) = 7.0d0
      call ensureflambda (-2000.0d0)
      call checkint ('ensureflambda low nflmda',
     &               nflmda,3005,nfail)
      call checkint ('ensureflambda low fli0',
     &               fli0,3003,nfail)
      call checkreal ('ensureflambda low shifted gkernel',
     &                gkernel(2,3003),7.0d0,1.0d-12,nfail)
c
c     low-side expansion also rebuilds osthead, ostnext and osthist
c
      call resetost (3,5,2)
      nosthist = 2
      call sethist (1,0.50d0,0.0d0,1.0d0,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,2.0d0,wlmda,wflmda)
      call buildostindex
      call ensureflambda (-2000.0d0)
c
c     after low-side expansion, fli0 shifts from 3 to 3003, so
c     flambda=0.0 now belongs to iflmda=3003
c
      ilmda = 2
      iflmda = 3003
      call ij_to_k (ilmda,iflmda,nlmda,k)
      call checkint ('ensureflambda rebuild osthist(1)',
     &               osthist(1),k,nfail)
      call checkint ('ensureflambda rebuild osthist(2)',
     &               osthist(2),k,nfail)
      call checkint ('ensureflambda rebuild osthead',
     &               osthead(ilmda,iflmda),2,nfail)
      call checkint ('ensureflambda rebuild ostnext head',
     &               ostnext(2),1,nfail)
      call checkint ('ensureflambda rebuild ostnext tail',
     &               ostnext(1),0,nfail)
      return
      end
c
c
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine testgkernels -- test g kernel routines        ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine testgkernels (nfail)
      use dlmda
      use math
      use ost
      implicit none
      integer nfail
      real*8 egbias,dgdl,dgdfl
      real*8 expected
      real*8 height
c
c
c     choose height so the normalized gaussian prefactor is one
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nosthist = 1
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
c
c     addgkernelhist from a single interior gaussian
c
      call addgkernelhist (1)
      call checkreal ('addgkernelhist center value',
     &                gkernel(3,3),1.0d0,1.0d-12,nfail)
      expected = exp(-1.0d0)
      call checkreal ('addgkernelhist off-center value',
     &                gkernel(4,4),expected,1.0d-12,nfail)
c
c     addgkernelhist includes left-boundary mirror image
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nosthist = 1
      call sethist (1,0.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      call addgkernelhist (1)
      call checkreal ('addgkernelhist left mirror center',
     &                gkernel(1,3),2.0d0,1.0d-12,nfail)
      expected = 2.0d0 * exp(-0.5d0)
      call checkreal ('addgkernelhist left mirror neighbor',
     &                gkernel(2,3),expected,1.0d-12,nfail)
c
c     addgkernelhist includes right-boundary mirror image
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nosthist = 1
      call sethist (1,1.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      call addgkernelhist (1)
      call checkreal ('addgkernelhist right mirror center',
     &                gkernel(5,3),2.0d0,1.0d-12,nfail)
      expected = 2.0d0 * exp(-0.5d0)
      call checkreal ('addgkernelhist right mirror neighbor',
     &                gkernel(4,3),expected,1.0d-12,nfail)
c
c     rebuild the original interior gaussian for later checks
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nosthist = 1
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
c
c     buildgkernel rebuilds the same full grid
c
      call buildgkernel
      call checkreal ('buildgkernel center value',
     &                gkernel(3,3),1.0d0,1.0d-12,nfail)
      expected = exp(-1.0d0)
      call checkreal ('buildgkernel off-center value',
     &                gkernel(4,4),expected,1.0d-12,nfail)
c
c     updategkernel adds only the newest gaussian to current grid
c
      gkernel = 0.0d0
      call updategkernel
      call checkreal ('updategkernel center value',
     &                gkernel(3,3),1.0d0,1.0d-12,nfail)
c
c     adding a taller gaussian updates only from the new history entry
c
      nosthist = 2
      call sethist (2,0.75d0,1.0d0,2.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
      expected = 2.0d0 + exp(-1.0d0)
      call checkreal ('updategkernel added new center',
     &                gkernel(4,4),expected,1.0d-12,nfail)
      expected = 1.0d0 + 2.0d0*exp(-1.0d0)
      call checkreal ('updategkernel keeps old plus new overlap',
     &                gkernel(3,3),expected,1.0d-12,nfail)
c
c     egkernel evaluates the continuous gaussian and derivatives
c
      ostlambda = 0.75d0
      ostdedl = 1.0d0
      call egkernel (egbias,dgdl,dgdfl)
      expected = 2.0d0 + exp(-1.0d0)
      call checkreal ('egkernel bias lambda=.75 flambda=1',
     &                egbias,expected,1.0d-12,nfail)
      call checkreal ('egkernel dgdl lambda=.75 flambda=1',
     &                dgdl,-4.0d0*exp(-1.0d0),1.0d-12,nfail)
      call checkreal ('egkernel dgdfl lambda=.75 flambda=1',
     &                dgdfl,-exp(-1.0d0),1.0d-12,nfail)
c
c     two gaussians in the same bin are both followed by ostnext
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nosthist = 2
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      ostlambda = 0.50d0
      ostdedl = 0.0d0
      call egkernel (egbias,dgdl,dgdfl)
      call checkreal ('egkernel same-bin linked bias',
     &                egbias,2.0d0,1.0d-12,nfail)
      call checkreal ('egkernel same-bin linked dgdl',
     &                dgdl,0.0d0,1.0d-12,nfail)
      call checkreal ('egkernel same-bin linked dgdfl',
     &                dgdfl,0.0d0,1.0d-12,nfail)
c
c     multiple bins are found through osthead lookup, including
c     two gaussians in one bin and one gaussian in another bin
c
      call resetost (5,5,4)
      height = 2.0d0 * pi * wlmda * wflmda
      nosthist = 3
      call sethist (1,0.50d0,0.0d0,1.0d0*height,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,2.0d0*height,wlmda,wflmda)
      call sethist (3,0.75d0,1.0d0,3.0d0*height,wlmda,wflmda)
      call buildostindex
      ostlambda = 0.75d0
      ostdedl = 1.0d0
      call egkernel (egbias,dgdl,dgdfl)
      expected = 3.0d0 + 3.0d0*exp(-1.0d0)
      call checkreal ('egkernel multi-bin mixed bias',
     &                egbias,expected,1.0d-12,nfail)
      call checkreal ('egkernel multi-bin mixed dgdl',
     &                dgdl,-12.0d0*exp(-1.0d0),1.0d-12,nfail)
      call checkreal ('egkernel multi-bin mixed dgdfl',
     &                dgdfl,-3.0d0*exp(-1.0d0),1.0d-12,nfail)
c
c     left endpoint includes both real and mirror gaussian images
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nosthist = 1
      call sethist (1,0.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      ostlambda = 0.0d0
      ostdedl = 0.0d0
      call egkernel (egbias,dgdl,dgdfl)
      call checkreal ('egkernel left mirror bias',
     &                egbias,2.0d0,1.0d-12,nfail)
      call checkreal ('egkernel left mirror dgdl',
     &                dgdl,0.0d0,1.0d-12,nfail)
c
c     right endpoint includes both real and mirror gaussian images
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nosthist = 1
      call sethist (1,1.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      ostlambda = 1.0d0
      ostdedl = 0.0d0
      call egkernel (egbias,dgdl,dgdfl)
      call checkreal ('egkernel right mirror bias',
     &                egbias,2.0d0,1.0d-12,nfail)
      call checkreal ('egkernel right mirror dgdl',
     &                dgdl,0.0d0,1.0d-12,nfail)
c
c     if dU/dlambda is outside the current flambda grid,
c     egkernel returns zero
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nosthist = 1
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      ostlambda = 0.50d0
      ostdedl = 100.0d0
      call egkernel (egbias,dgdl,dgdfl)
      call checkreal ('egkernel outside flambda bias',
     &                egbias,0.0d0,1.0d-12,nfail)
      call checkreal ('egkernel outside flambda dgdl',
     &                dgdl,0.0d0,1.0d-12,nfail)
      call checkreal ('egkernel outside flambda dgdfl',
     &                dgdfl,0.0d0,1.0d-12,nfail)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine testfkernel -- test f kernel construction     ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine testfkernel (nfail)
      use bath
      use math
      use ost
      use units
      implicit none
      integer nfail
      real*8 rt
      real*8 expected
      real*8 height
      real*8 w1,w2,w3,w4,w5
c
c     two nonzero gkernel weights at flambda=-1 and +1
c
      call resetost (5,5,1)
      rt = gasconst * kelvin
      gkernel(3,2) = log(2.0d0) * rt
      gkernel(3,4) = log(4.0d0) * rt
      call buildfkernel
      expected = 1.0d0 / 3.0d0
      call checkreal ('buildfkernel weighted mean force',
     &                fkernel(3),expected,1.0d-12,nfail)
      call checkreal ('buildfkernel empty row',
     &                fkernel(1),0.0d0,1.0d-12,nfail)
c
c     build gkernel incrementally from multiple gaussians using
c     updategkernel, then construct the f kernel
c
      call resetost (5,5,4)
      rt = gasconst * kelvin
      height = 2.0d0 * pi * wlmda * wflmda
c
      nosthist = 1
      call sethist (1,0.25d0,-1.0d0,1.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
c
      nosthist = 2
      call sethist (2,0.50d0, 0.0d0,2.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
c
      nosthist = 3
      call sethist (3,0.75d0, 1.0d0,3.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
c
      call buildfkernel
c
c     expected weighted mean at lambda=0.5
c
      w1 = exp(gkernel(3,1)/rt)
      w2 = exp(gkernel(3,2)/rt)
      w3 = exp(gkernel(3,3)/rt)
      w4 = exp(gkernel(3,4)/rt)
      w5 = exp(gkernel(3,5)/rt)
      
      expected = (-2.d0*w1 - w2 + w4 + 2.d0*w5)
     &     /(w1+w2+w3+w4+w5)
      call checkreal ('buildfkernel from updated gkernel',
     &                fkernel(3),expected,1.0d-12,nfail)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine testefkernel -- test DeltaG interpolation     ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine testefkernel (nfail)
      use ost
      implicit none
      integer nfail
      integer i
      real*8 eostlmda,dfdl
      real*8 expected
c
c
c     use fkernel(lambda)=lambda and integrate to lambda=0.375
c
      call resetost (5,5,1)
      do i = 1, nlmda
         fkernel(i) = dble(i-1) * wlmda
      end do
      ostlambda = 0.375d0
      call efkernel (eostlmda,dfdl)
      expected = 0.5d0 * ostlambda * ostlambda
      call checkreal ('efkernel DeltaG lambda=.375',
     &                eostlmda,expected,1.0d-12,nfail)
      call checkreal ('efkernel dDeltaG/dlambda',
     &                dfdl,ostlambda,1.0d-12,nfail)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine testmapsub -- test lambda mapping schemes     ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine testmapsub (nfail)
      use dlmda
      use mutant
      use ost
      implicit none
      integer nfail
c
c
c     test exponential sublambda maps and chain rule derivatives
c
      call resetost (5,5,1)
      ostlambda = 0.25d0
      ostpmap = 'EXP'
      ostemap = 'EXP'
      ostvmap = 'EXP'
      ostepexp = 2
      ostemexp = 3
      ostevexp = 4
      call mapsublmda
      call checkreal ('mapsublmda exponential plambda',
     &                plambda,0.0625d0,1.0d-12,nfail)
      call checkreal ('mapsublmda exponential dpldlmda',
     &                dpldlmda,0.5d0,1.0d-12,nfail)
      call checkreal ('mapsublmda exponential d2pldlmda2',
     &                d2pldlmda2,2.0d0,1.0d-12,nfail)
      call checkreal ('mapsublmda exponential elambda',
     &                elambda,0.015625d0,1.0d-12,nfail)
      call checkreal ('mapsublmda exponential deldlmda',
     &                deldlmda,0.1875d0,1.0d-12,nfail)
      call checkreal ('mapsublmda exponential d2eldlmda2',
     &                d2eldlmda2,1.5d0,1.0d-12,nfail)
      call checkreal ('mapsublmda exponential vlambda',
     &                vlambda,0.00390625d0,1.0d-12,nfail)
      call checkreal ('mapsublmda exponential dvldlmda',
     &                dvldlmda,0.0625d0,1.0d-12,nfail)
      call checkreal ('mapsublmda exponential d2vldlmda2',
     &                d2vldlmda2,0.75d0,1.0d-12,nfail)
c
c     test shifted inverse-power sublambda maps and derivatives
c
      call resetost (5,5,1)
      ostlambda = 0.25d0
      ostpmap = 'INV'
      ostemap = 'INV'
      ostvmap = 'INV'
      ostinvepn = 2
      ostinvemn = 3
      ostinvevn = 4
      ostinvepeps = 0.01d0
      ostinvemeps = 0.02d0
      ostinveveps = 0.03d0
      call mapsublmda
      call checkreal ('mapsublmda invpower plambda',
     &                plambda,0.452936557937477d0,1.0d-12,nfail)
      call checkreal ('mapsublmda invpower dpldlmda',
     &                dpldlmda,1.08352945028593d0,1.0d-12,nfail)
      call checkreal ('mapsublmda invpower d2pldlmda2',
     &                d2pldlmda2,-2.08371048131911d0,1.0d-12,
     &                nfail)
      call checkreal ('mapsublmda invpower elambda',
     &                elambda,0.509927040983045d0,1.0d-12,nfail)
      call checkreal ('mapsublmda invpower deldlmda',
     &                deldlmda,1.08536378202464d0,1.0d-12,nfail)
      call checkreal ('mapsublmda invpower d2eldlmda2',
     &                d2eldlmda2,-2.67991057290036d0,1.0d-12,
     &                nfail)
      call checkreal ('mapsublmda invpower vlambda',
     &                vlambda,0.526434441030379d0,1.0d-12,nfail)
      call checkreal ('mapsublmda invpower dvldlmda',
     &                dvldlmda,1.09852311505228d0,1.0d-12,nfail)
      call checkreal ('mapsublmda invpower d2vldlmda2',
     &                d2vldlmda2,-2.94247262960432d0,1.0d-12,
     &                nfail)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine testmeta -- test metadynamics bias helpers    ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine testmeta (nfail)
      use math
      use ost
      implicit none
      integer nfail
      real*8 vbias,dvdl
      real*8 pref
      real*8 expected
      real*8 metadeltag
c
c
c     one normalized 1D gaussian centered at lambda=0.5
c
      call resetost (5,5,1)
      call resetmeta (2)
      nmetahist = 1
      metalhist(1) = 0.5d0
      metahhist(1) = 2.0d0
      metawhist(1) = 0.25d0
      pref = metahhist(1) / (metawhist(1)*sqrt(2.0d0*pi))
      call emetabias (0.5d0,vbias,dvdl)
      call checkreal ('emetabias center value',
     &                vbias,pref,1.0d-12,nfail)
      call checkreal ('emetabias center derivative',
     &                dvdl,0.0d0,1.0d-12,nfail)
      call emetabias (0.75d0,vbias,dvdl)
      expected = pref * exp(-0.5d0)
      call checkreal ('emetabias off-center value',
     &                vbias,expected,1.0d-12,nfail)
      call checkreal ('emetabias off-center derivative',
     &                dvdl,-4.0d0*expected,1.0d-12,nfail)
c
c     symmetric gaussian has zero endpoint free energy difference
c
      call checkreal ('metadeltag symmetric gaussian',
     &                metadeltag(),0.0d0,1.0d-12,nfail)
c
c     resizing preserves old metadynamics history and zeros new slots
c
      nmetahist = 2
      metalhist(2) = 0.25d0
      metahhist(2) = 3.0d0
      metawhist(2) = 0.125d0
      call resizemeta
      call checkint ('resizemeta size input size=2',
     &               sizemetahist,4,nfail)
      call checkreal ('resizemeta preserve metalhist(2)',
     &                metalhist(2),0.25d0,1.0d-12,nfail)
      call checkreal ('resizemeta init metalhist(3)',
     &                metalhist(3),0.0d0,1.0d-12,nfail)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine resetost -- allocate and initialize ost state  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine resetost (nl,nf,nhist)
      use ost
      implicit none
      integer nl,nf,nhist
      integer i,j
c
c
c     clear any previous allocation
c
      if (allocated(osthhist))  deallocate (osthhist)
      if (allocated(osthist))  deallocate (osthist)
      if (allocated(osthead))  deallocate (osthead)
      if (allocated(ostnext))  deallocate (ostnext)
      if (allocated(ostlhist))  deallocate (ostlhist)
      if (allocated(ostfhist))  deallocate (ostfhist)
      if (allocated(ostwlhist))  deallocate (ostwlhist)
      if (allocated(ostwfhist))  deallocate (ostwfhist)
      if (allocated(fkernel))  deallocate (fkernel)
      if (allocated(gkernel))  deallocate (gkernel)
      if (allocated(metalhist))  deallocate (metalhist)
      if (allocated(metahhist))  deallocate (metahhist)
      if (allocated(metawhist))  deallocate (metawhist)
c
c     set scalar state
c
      nlmda = nl
      nflmda = nf
      fli0 = (nflmda + 1) / 2
      wlmda = 1.0d0 / dble(nlmda-1)
      wflmda = 1.0d0
      wlmda2 = 0.5d0 * wlmda
      wflmda2 = 0.5d0 * wflmda
      nosthist = 0
      sizeosthist = nhist
      ostlambda = 0.0d0
      ostlambdaavg = 0.0d0
      ostdedl = 0.0d0
      ostdedlavg = 0.0d0
      deffdl = 0.0d0
      ostpmap = 'QNT'
      ostemap = 'QNT'
      ostvmap = 'QNT'
      ostepexp = 1
      ostemexp = 1
      ostevexp = 1
      ostinvepn = 1
      ostinvemn = 1
      ostinvevn = 1
      ostinvepeps = 0.0d0
      ostinvemeps = 0.0d0
      ostinveveps = 0.0d0
      hbias = 0.0d0
      eosttot = 0.0d0
      oststdev = 1.0d0
c
c     allocate arrays
c
      allocate (osthist(sizeosthist))
      allocate (osthead(nlmda,nflmda))
      allocate (ostnext(sizeosthist))
      allocate (ostlhist(sizeosthist))
      allocate (ostfhist(sizeosthist))
      allocate (osthhist(sizeosthist))
      allocate (ostwlhist(sizeosthist))
      allocate (ostwfhist(sizeosthist))
      allocate (fkernel(nlmda))
      allocate (gkernel(nlmda,nflmda))
c
c     initialize arrays
c
      do i = 1, sizeosthist
         osthist(i) = 0
         ostnext(i) = 0
         ostlhist(i) = 0.0d0
         ostfhist(i) = 0.0d0
         osthhist(i) = 0.0d0
         ostwlhist(i) = 0.0d0
         ostwfhist(i) = 0.0d0
      end do
      do i = 1, nlmda
         fkernel(i) = 0.0d0
         do j = 1, nflmda
            osthead(i,j) = 0
            gkernel(i,j) = 0.0d0
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine resetmeta -- allocate metadynamics history    ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine resetmeta (nhist)
      use ost
      implicit none
      integer nhist
      integer i
c
c
c     clear and allocate metadynamics history arrays
c
      if (allocated(metalhist))  deallocate (metalhist)
      if (allocated(metahhist))  deallocate (metahhist)
      if (allocated(metawhist))  deallocate (metawhist)
      sizemetahist = nhist
      nmetahist = 0
      allocate (metalhist(sizemetahist))
      allocate (metahhist(sizemetahist))
      allocate (metawhist(sizemetahist))
      do i = 1, sizemetahist
         metalhist(i) = 0.0d0
         metahhist(i) = 0.0d0
         metawhist(i) = 0.0d0
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine sethist -- set one gaussian history entry     ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine sethist (ihist,lambda,flmda,height,sigl,sigf)
      use ost
      implicit none
      integer ihist
      integer ilmda,iflmda,k
      integer lambdabin,flambdabin
      real*8 lambda,flmda
      real*8 height,sigl,sigf
c
c
c     save real gaussian center, parameters and packed lookup bin
c
      ilmda = lambdabin(lambda)
      iflmda = flambdabin(flmda)
      call ij_to_k (ilmda,iflmda,nlmda,k)
      osthist(ihist) = k
      ostlhist(ihist) = lambda
      ostfhist(ihist) = flmda
      osthhist(ihist) = height
      ostwlhist(ihist) = sigl
      ostwfhist(ihist) = sigf
      ostnext(ihist) = 0
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine checkint -- check integer expected output     ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine checkint (name,actual,expected,nfail)
      implicit none
      integer actual,expected,nfail
      character*(*) name
c
c
      if (actual .eq. expected) then
         write (*,10)  name,expected,actual
   10    format (' PASS ',a,': expected ',i12,' actual ',i12)
      else
         nfail = nfail + 1
         write (*,20)  name,expected,actual
   20    format (' FAIL ',a,': expected ',i12,' actual ',i12)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine checkreal -- check real expected output       ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine checkreal (name,actual,expected,tol,nfail)
      implicit none
      integer nfail
      real*8 actual,expected,tol
      real*8 diff
      character*(*) name
c
c
      diff = abs(actual-expected)
      if (diff .le. tol) then
         write (*,10)  name,expected,actual
   10    format (' PASS ',a,': expected ',g24.15,' actual ',g24.15)
      else
         nfail = nfail + 1
         write (*,20)  name,expected,actual,diff
   20    format (' FAIL ',a,': expected ',g24.15,' actual ',
     &           g24.15,' diff ',g12.5)
      end if
      return
      end
