c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtpdb  --  output of Protein Data Bank file  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtpdb" writes out a set of Protein Data Bank coordinates
c     to an external disk file
c
c
      subroutine prtpdb (ipdb)
      implicit none
      include 'sizes.i'
      include 'files.i'
      include 'pdb.i'
      include 'sequen.i'
      include 'titles.i'
      integer i,k,ipdb
      integer start,stop,resnumb
      integer resid(maxres)
      logical opened
      character*1 chnname
      character*1 chain(maxres)
      character*3 resname
      character*120 pdbfile
c
c
c     open output unit if not already done
c
      inquire (unit=ipdb,opened=opened)
      if (.not. opened) then
         pdbfile = filename(1:leng)//'.pdb'
         call version (pdbfile,'new')
         open (unit=ipdb,file=pdbfile,status='new')
      end if
c
c     write out the header lines and the title
c
      if (ltitle .eq. 0) then
         write (ipdb,10)
   10    format ('HEADER',/,'COMPND',/,'SOURCE')
      else
         write (ipdb,20)  title(1:ltitle)
   20    format ('HEADER',4x,a,/,'COMPND',/,'SOURCE')
      end if
c
c     find the chain name and chain position for each residue
c
      do i = 1, nchain
         start = ichain(1,i)
         stop = ichain(2,i)
         do k = start, stop
            resid(k) = k - start + 1
            chain(k) = chnnam(i)
         end do
      end do
c
c     change some TINKER residue names to match PDB standards
c
      do i = 1, npdb
         if (resnam(i) .eq. 'CYX')  resnam(i) = 'CYS'
         if (resnam(i) .eq. 'HID')  resnam(i) = 'HIS'
         if (resnam(i) .eq. 'HIE')  resnam(i) = 'HIS'
         if (resnam(i) .eq. 'HIP')  resnam(i) = 'HIS'
      end do
c
c     next, write the coordinates for each PDB atom
c
      do i = 1, npdb
         resname = resnam(i)
         if (resname(2:3) .eq. '  ')  resname = '  '//resname(1:1)
         if (resname(3:3) .eq. ' ')  resname = ' '//resname(1:2)
         if (pdbtyp(i) .eq. 'ATOM  ') then
            resnumb = resid(resnum(i))
            chnname = chain(resnum(i))
         else
            resnumb = resnum(i)
            chnname = ' '
         end if
         write (ipdb,30)  pdbtyp(i),i,atmnam(i),resname,chnname,
     &                    resnumb,xpdb(i),ypdb(i),zpdb(i)
   30    format (a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)
      end do
c
c     finally, write any connectivity records for PDB atoms
c
      do i = 1, npdb
         if (npdb12(i) .ne. 0) then
            write (ipdb,40)  i,(ipdb12(k,i),k=1,npdb12(i))
   40       format ('CONECT',5i5)
         end if
      end do
      write (ipdb,50)
   50 format ('END')
c     close (unit=ipdb)
      return
      end
