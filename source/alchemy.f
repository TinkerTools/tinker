c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 1991 by Shawn Huston & Jay William Ponder  ##
c     ##                    All Rights Reserved                    ##
c     ###############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program alchemy  --  perform free energy perturbation  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "alchemy" computes the free energy difference corresponding
c     to a small perturbation by Boltzmann weighting the potential
c     energy difference over a number of sample states; current
c     version (incorrectly) considers the charge energy to be
c     intermolecular in finding the perturbation energies
c
c     variables and parameters :
c
c     nlamb   number of lambda values for energy computation
c     delta   step size for the perturbation in lambda
c     nrun    number of steps over which to calculate averages
c
c     deplus       energy change for postive delta lambda
c     deminus      energy change for negative delta lambda
c     sdep,sdem    accumulated energy changes
c     adep,adem    (running) average energy changes
c     adapb,adamb  average over block averages of free energy changes
c     bdep,bdem    (block) accumulated energy changes
c     badep,badem  (block) average energy changes
c     vdep,vdem    SD of energy changes from variance of subaverages
c     fdep,fdem    fluctuations of perturbation energies
c     se_ap,se_am  standard error of free energy changes
c
c
      program alchemy
      implicit none
      include 'sizes.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'energi.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'katoms.i'
      include 'mutant.i'
      include 'potent.i'
      include 'units.i'
      include 'usage.i'
      integer maxblock,maxstep
      parameter (maxblock=100)
      parameter (maxstep=5000)
      integer i,j,k,ixyz,start,stop
      integer nblock,nrun
      integer modrun,modblock
      integer istep,nstep,ilamb,nlamb
      integer lext,next,freeunit
      real*8 delta,lam0,lamp,lamm
      real*8 rt,e0,energy,pos,neg,temp
      real*8 eplus,eminus
      real*8 deplus,deminus
      real*8 sdep,s2dep,sdep2
      real*8 sdem,s2dem,sdem2
      real*8 spos,bpos,sdapb
      real*8 sneg,bneg,sdamb
      real*8 bdep,adep,a2dep,adep2
      real*8 bdem,adem,a2dem,adem2
      real*8 da,dap,dam,adapb,adamb
      real*8 v,vdap,vdam,vdep,vdem
      real*8 fdep,fdem,bda
      real*8 se_ep,se_em,se_ap,se_am
      real*8 eb0,ebpos,ebneg
      real*8 ea0,eapos,eaneg
      real*8 eit0,eitpos,eitneg
      real*8 et0,etpos,etneg
      real*8 ev0,evpos,evneg
      real*8 ec0,ecpos,ecneg
      real*8 sbp,sap,sitp,stp,svp,scp
      real*8 sbm,sam,sitm,stm,svm,scm
      real*8 abp,aap,aitp,atp,avp,acp
      real*8 abm,aam,aitm,atm,avm,acm
      real*8 badep(maxblock),badem(maxblock)
      real*8 bdap(maxblock),bdam(maxblock)
      real*8 nrg(3,maxstep),cb(3,maxstep)
      real*8 ca(3,maxstep),cit(3,maxstep)
      real*8 ct(3,maxstep),cv(3,maxstep)
      real*8 cc(3,maxstep)
      logical dogeom
      character*1 answer
      character*7 ext
      character*120 xyzfile
      character*120 record
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     get the numbers of the files to be used
c
      start = 0
      stop = 0
      write (iout,10)
   10 format (/,' Numbers of First and Last File to Analyze :  ',$)
      read (input,20)  record
   20 format (a120)
      read (record,*,err=30,end=30)  start,stop
   30 continue
      if (start .eq. 0)  start = 1
      if (stop .eq. 0)  stop = start
      nrun  = stop - start + 1
c
c     obtain the lambda values to be calculated
c
      delta = 0.0d0
      write (iout,40)
   40 format (/,' Enter the Lambda Increment for FEP :  ',$)
      read (input,50)  delta
   50 format (f7.4)
      nlamb = 3
      lam0 = lambda
      lamp = min(1.0d0,lambda+delta)
      lamm = max(0.0d0,lambda-delta)
c
c     obtain the target temperature value
c
      temp = 0.0d0
      write (iout,60)
   60 format (/,' Enter the System Temperature [300 K] :  ',$)
      read (input,70)  temp
   70 format (f20.0)
      if (temp .eq. 0.0d0)  temp = 300.0d0
      rt = gasconst * temp
c
c     set number of steps for running averages and block averages
c
      nblock = 0
      write (iout,80)
   80 format (/,' Enter Number of Blocks for Sub-Averages [1] :  ',$)
      read (input,90)  nblock
   90 format (i10)
      if (nblock .eq. 0)  nblock = 1
      nblock = nrun / nblock
c
c     decide whether to include the intramolecular energies
c
      dogeom = .true.
      write (iout,100)
  100 format (/,' Consider only Intermolecular Perturbation',
     &           ' Energy [N] :  ',$)
      read (input,110)  record
  110 format (a120)
      next = 1
      call gettext (record,answer,next)
      call upcase (answer)
      if (answer .eq. 'Y')  dogeom = .false.
      if (dogeom) then
         write (iout,120)
  120    format (/,' Calculation will Involve Full Perturbation',
     &              ' Energy ')
      else
         write (iout,130)
  130    format (/,' Calculation will Consider Only Intermolecular',
     &              ' Interactions ')
      end if
c
c     cycle over the coordinate files once per lambda value
c
      do ilamb = 1, nlamb
         i = start
         istep = 0
         if (ilamb .eq. 2)  lambda = lamp
         if (ilamb .eq. 3)  lambda = lamm
         call hybrid
c
c     read in the next molecular dynamics coordinate frame
c
         do while (i.ge.start .and. i.le.stop)
            istep = istep + 1
            lext = 3
            call numeral (i,ext,lext)
            xyzfile = filename(1:leng)//'.'//ext(1:lext)
            call version (xyzfile,'old')
            ixyz = freeunit ()
            open (unit=ixyz,file=xyzfile,status='old',err=160)
            call readxyz (ixyz)
            close (unit=ixyz)
            call hatom
c
c     select interactions for perturbation energy calculation
c
            do j = 1, n
               use(j) = .false.
            end do
            do j = 1, nmut
               use(imut(j)) = .true.
            end do
            if (.not. dogeom) then
               use_bond = .false.
               use_angle = .false.
               use_strbnd = .false.
               use_urey = .false.
               use_imptor = .false.
               use_tors = .false.
               use_strtor = .false.
            end if
c
c     compute and store energy components for the current lambda
c
            nrg(ilamb,istep) = energy ()
            cb(ilamb,istep) = eb
            ca(ilamb,istep) = ea
            cit(ilamb,istep) = eit
            ct(ilamb,istep) = et
            cv(ilamb,istep) = ev
            cc(ilamb,istep) = ec
            if (verbose) then
               if (istep .eq. 1) then
                  write (iout,140)
  140             format (/,2x,'Step',6x,'EB',8x,'EA',7x,'EIT',
     &                       8x,'ET',8x,'EV',8x,'EC',/)
               end if
               write (iout,150)  istep,eb,ea,eit,et,ev,ec
  150          format (i6,6f10.4)
            end if
  160       continue
            i = i + 1
         end do
      end do
      nstep = istep
c
c     get free energy change by averaging over all frames
c
      do istep = 1, nstep
         e0 = nrg(1,istep)
         eplus = nrg(2,istep)
         eminus = nrg(3,istep)
         ev0 = cv(1,istep)
         evpos = cv(2,istep)
         evneg = cv(3,istep)
         ec0 = cc(1,istep)
         ecpos = cc(2,istep)
         ecneg = cc(3,istep)
         if (dogeom) then
            eb0 = cb(1,istep)
            ebpos = cb(2,istep)
            ebneg = cb(3,istep)
            ea0 = ca(1,istep)
            eapos = ca(2,istep)
            eaneg = ca(3,istep)
            eit0 = cit(1,istep)
            eitpos = cit(2,istep)
            eitneg = cit(3,istep)
            et0 = ct(1,istep)
            etpos = ct(2,istep)
            etneg = ct(3,istep)
         end if
         modrun = mod(istep-1,nrun)
         modblock = mod(istep-1,nblock)
c
c     zero out summation variables for new running average
c
         if (modrun .eq. 0) then
            sdep = 0.0d0
            s2dep = 0.0d0
            sdep2 = 0.0d0
            sdem = 0.0d0
            s2dem = 0.0d0
            sdem2 = 0.0d0
            spos = 0.0d0
            sneg = 0.0d0
            sdapb = 0.0d0
            sdamb = 0.0d0
            sbp = 0.0d0
            sbm = 0.0d0
            sap = 0.0d0
            sam = 0.0d0
            sitp = 0.0d0
            sitm = 0.0d0
            stp = 0.0d0
            stm = 0.0d0
            svp = 0.0d0
            svm = 0.0d0
            scp = 0.0d0
            scm = 0.0d0
            fdep = 0.0d0
            fdem = 0.0d0
            vdep = 0.0d0
            vdem = 0.0d0
            vdap = 0.0d0
            vdam = 0.0d0
            k = 0
         end if
c
c     zero out summation variables for new block average
c
         if (modblock .eq. 0) then
            bdep = 0.0d0
            bdem = 0.0d0
            bpos = 0.0d0
            bneg = 0.0d0
         end if
         modrun = mod(istep,nrun)
         modblock = mod(istep,nblock)
c
c     accumulate statistics
c
         deplus = eplus - e0
         deminus = eminus - e0
         if (verbose) then
            if (modblock.eq.1 .or. nblock.eq.1) then
               write (iout,170)
  170          format (/,2x,'Step',8x,'E0',10x,'EP',10x,'EM',
     &                    10x,'DEP',9x,'DEM',/)
            end if
            write (iout,180)  istep,e0,eplus,eminus,deplus,deminus
  180       format (i6,5f12.4)
         end if
         pos = exp(-deplus/rt)
         neg = exp(-deminus/rt)
         bdep = bdep + deplus
         bdem = bdem + deminus
         bpos = bpos + pos
         bneg = bneg + neg
         sdep = sdep + deplus
         sdem = sdem + deminus
         s2dep = s2dep + deplus*deplus
         s2dem = s2dem + deminus*deminus
         spos = spos + pos
         sneg = sneg + neg
         svp = svp + evpos - ev0
         svm = svm + evneg - ev0
         scp = scp + ecpos - ec0
         scm = scm + ecneg - ec0
         if (dogeom) then
            sbp = sbp + ebpos - eb0
            sbm = sbm + ebneg - eb0
            sap = sap + eapos - ea0
            sam = sam + eaneg - ea0
            sitp = sitp + eitpos - eit0
            sitm = sitm + eitneg - eit0
            stp = stp + etpos - et0
            stm = stm + etneg - et0
         end if
c
c     calculate block averages
c
         if (modblock .eq. 0) then
            k = k + 1
            badep(k) = bdep / dble(nblock)
            badem(k) = bdem / dble(nblock)
            bda = bpos / dble(nblock)
            bda = -rt * log(bda)
            bdap(k) = bdap(k) +  bda
            bda = bneg / dble(nblock)
            bda = -rt * log(bda)
            bdam(k) = bdam(k) + bda
            sdapb = sdapb + bdap(k)
            sdamb = sdamb + bdam(k)
            if (verbose .or. k.eq.1) then
               write (iout,190)
  190          format (/,2x,'Block',6x,'NStep',7x,'BADEP',
     &                    7x,'BADEM',8x,'BDAP',8x,'BDAM',/)
            end if
            write (iout,200)  k,istep,badep(k),badem(k),
     &                        bdap(k),bdam(k)
  200       format (i6,5x,i6,1x,4f12.4)
         end if
c
c     calculate running averages for potential energy
c
         if (modrun .eq. 0) then
            adep = sdep / dble(nrun)
            adem = sdem / dble(nrun)
            a2dep = s2dep / dble(nrun)
            a2dem = s2dem / dble(nrun)
            adep2 = adep * adep
            adem2 = adem * adem
            fdep =  sqrt(a2dep - adep2)
            fdem =  sqrt(a2dem - adem2)
            do k = 1, nrun / nblock
               v = (badep(k) - adep)**2
               vdep = vdep + v
               v = (badem(k) - adem)**2
               vdem = vdem + v
            end do
            vdep = vdep / dble(nrun/nblock)
            se_ep = sqrt(vdep / dble(nrun/nblock))
            vdem = vdem / dble(nrun/nblock)
            se_em = sqrt(vdem / dble(nrun/nblock))
c
c     calculate running averages for free energy
c
            da = spos / dble(nrun)
            da = -rt * log (da)
            dap = da
            da = sneg / dble(nrun)
            da = -rt * log (da)
            dam = da
            adapb = sdapb / dble(nrun/nblock)
            adamb = sdamb / dble(nrun/nblock)
            do k = 1, nrun/nblock
               v = (bdap(k) - adapb)**2
               vdap = vdap + v
               v = (bdam(k) - adamb)**2
               vdam = vdam + v
            end do
            vdap = vdap / dble(nrun/nblock)
            se_ap = sqrt(vdap / dble(nrun/nblock))
            vdam = vdam / dble(nrun/nblock)
            se_am = sqrt(vdam / dble(nrun/nblock))
c
c     calculate running averages for energy components
c
            avp = svp / dble(nrun)
            avm = svm / dble(nrun)
            acp = scp / dble(nrun)
            acm = scm / dble(nrun)
            if (dogeom) then
               abp = sbp / dble(nrun)
               abm = sbm / dble(nrun)
               aap = sap / dble(nrun)
               aam = sam / dble(nrun)
               aitp = sitp / dble(nrun)
               aitm = sitm / dble(nrun)
               atp = stp / dble(nrun)
               atm = stm / dble(nrun)
            end if
            sdep = 0.0d0
            sdem = 0.0d0
c
c     write information about running averages and block averages
c
            write (iout,210)  nstep,nrun/nblock
  210       format (/,' Running Averages over',i5,' Steps',
     &                 ' with Std Error from',i4,' Blocks :')
            write (iout,220)
  220       format (/,' Free Energy :')
            write (iout,230)  dap,se_ap
  230       format (/,' DA(+) =',f12.4,' with Std Error',f10.4)
            write (iout,240)  dam,se_am
  240       format (' DA(-) =',f12.4,' with Std Error',f10.4)
            write (iout,250)
  250       format (/,' Potential Energy :')
            write (iout,260)  adep,fdep,se_ep
  260       format (/,' DE(+) =',f12.4,' with Fluct',f10.4,
     &                 ' and Std Error',f10.4)
            write (iout,270)  adem,fdem,se_em
  270       format (' DE(-) =',f12.4,' with Fluct',f10.4,
     &                 ' and Std Error',f10.4)
            write (iout,280)
  280       format (/,' Component Energies :',/)
            if (dogeom) then
               write (iout,290)  abp,abm
  290          format (' BOND  +/- :',f12.4,5x,f12.4)
               write (iout,300)  aap,aam
  300          format (' ANGLE +/- :',f12.4,5x,f12.4)
               write (iout,310)  aitp,aitm
  310          format (' IMPT  +/- :',f12.4,5x,f12.4)
               write (iout,320)  atp,atm
  320          format (' TORS  +/- :',f12.4,5x,f12.4)
            end if
            write (iout,330)  avp,avm
  330       format (' VDW   +/- :',f12.4,5x,f12.4)
            write (iout,340)  acp,acm
  340       format (' CHG   +/- :',f12.4,5x,f12.4)
         end if
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
