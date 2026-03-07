c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2025  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine readcif  --  input of PDBx/mmCIF format file  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readcif" gets a set of coordinates in RCSB PDBx/mmCIF
c     format from an external file
c
c
      subroutine readcif (icif)
      use boxes
      use files
      use inform
      use iounit
      use pdb
      use resdue
      use sequen
      use titles
      implicit none
      integer i,j,k,icif
      integer start,stop
      integer index,serial
      integer next,nxtlast
      integer residue,reslast
      integer model
      integer trimtext
      real*8 xx,yy,zz
      real*8 xbx,ybx,zbx
      real*8 aan,ban,gan
      real*8 occupy,bfac
      logical exist,opened
      logical first
      character*1 letter
      character*1 chain,chnlast
      character*1 altloc,formal
      character*1 insert,inslast
      character*1, allocatable :: chnatm(:)
      character*3 resname,atmsymb
      character*3 namelast
      character*4 atmname
      character*6 remark
      character*20 float
      character*240 ciffile
      character*240 record
      character*240 string
      save first
      data first  / .true. /
c
c
c     open the input file if it has not already been done
c
      inquire (unit=icif,opened=opened)
      if (.not. opened) then
         ciffile = filename(1:leng)//'.cif'
         call version (ciffile,'old')
         inquire (file=ciffile,exist=exist)
         if (exist) then
            open (unit=icif,file=ciffile,status='old')
            rewind (unit=icif)
         else
            write (iout,10)
   10       format (/,' READCIF  --  Unable to Find the PDBx/mmCIF',
     &                 ' Format File')
            call fatal
         end if
      end if
c
c     get alternate sites, chains and insertions to be used
c
      if (first)  call scancif (icif)
c
c     initialize title, atom and residue counters and name
c
      title = ' '
      ltitle = 0
      npdb = 0
      nres = 0
      reslast = maxres
      namelast = '   '
      chnlast = ' '
c
c     perform dynamic allocation of some local arrays
c
      allocate (chnatm(maxatm))
c
c     perform dynamic allocation of some global arrays
c
      if (first) then
         first = .false.
         if (.not. allocated(resnum))  allocate (resnum(maxatm))
         if (.not. allocated(resatm))  allocate (resatm(2,maxatm))
         if (.not. allocated(pdbmod))  allocate (pdbmod(maxatm))
         if (.not. allocated(xpdb))  allocate (xpdb(maxatm))
         if (.not. allocated(ypdb))  allocate (ypdb(maxatm))
         if (.not. allocated(zpdb))  allocate (zpdb(maxatm))
         if (.not. allocated(pdbres))  allocate (pdbres(maxatm))
         if (.not. allocated(pdbsym))  allocate (pdbsym(maxatm))
         if (.not. allocated(pdbatm))  allocate (pdbatm(maxatm))
         if (.not. allocated(pdbrec))  allocate (pdbrec(maxatm))
      end if
c
c     process info and individual atoms from PDBx/mmCIF file
c
      do while (.true.)
         read (icif,20,err=110,end=110)  record
   20    format (a240)
         call upcase (record)
         remark = record(1:6)
         if (record(1:13) .eq. '_STRUCT.TITLE') then
            next = 14
            call getstring (record,title,next)
            ltitle = trimtext (title)
            if (ltitle .eq. 0) then
               read (icif,30,err=40,end=40)  record
   30          format (a240)
               call upcase (record)
               next = 1
               call getstring (record,title,next)
               ltitle = trimtext (title)
            end if
   40       continue
         else if (record(1:21) .eq. '_DATABASE_2.PDBX_DOI') then
            if (ltitle .eq. 0) then
               next = 21
               ltitle = trimtext (record)
               title = record(21:ltitle)
               ltitle = ltitle - 20
               if (ltitle .eq. 0) then
                  read (icif,50,err=60,end=60)  record
   50             format (a240)
                  call upcase (record)
                  ltitle = trimtext (record)
                  title = record(1:ltitle)
               end if
   60          continue
            end if
         else if (record(1:14) .eq. '_CELL.LENGTH_A') then
            next = 15
            call getfloat (record,xbx,next)
         else if (record(1:14) .eq. '_CELL.LENGTH_B') then
            next = 15
            call getfloat (record,ybx,next)
         else if (record(1:14) .eq. '_CELL.LENGTH_C') then
            next = 15
            call getfloat (record,zbx,next)
         else if (record(1:17) .eq. '_CELL.ANGLE_ALPHA') then
            next = 18
            call getfloat (record,aan,next)
         else if (record(1:16) .eq. '_CELL.ANGLE_BETA') then
            next = 17
            call getfloat (record,ban,next)
         else if (record(1:17) .eq. '_CELL.ANGLE_GAMMA') then
            next = 18
            call getfloat (record,gan,next)
         else if (record(1:5) .eq. 'ATOM ') then
            remark = 'ATOM  '
            next = 6
            call getnumb (record,serial,next)
            call gettext (record,atmsymb,next)
            call gettext (record,atmname,next)
            call gettext (record,altloc,next)
            call gettext (record,resname,next)
            call gettext (record,chain,next)
            call gettext (record,letter,next)
            call getnumb (record,residue,next)
            call gettext (record,insert,next)
            call getfloat (record,xx,next)
            call getfloat (record,yy,next)
            call getfloat (record,zz,next)
            call getfloat (record,occupy,next)
            call getfloat (record,bfac,next)
            call gettext (record,formal,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call getnumb (record,model,next)
            if (altloc .eq. '.')  altloc = ' '
            if (insert .eq. '?')  insert = ' '
            if (formal .eq. '?')  formal = ' '
            if (index(chnsym,chain) .eq. 0)  goto 80
            if (altloc.ne.' ' .and. altloc.ne.altsym)  goto 80
            if (insert.ne.' ' .and. index(instyp,insert).eq.0)  goto 80
            if (model .ne. imodel)  goto 80
            call fixcif (resname,atmname)
            if (resname .eq. 'HOH') then
               remark = 'HETATM'
            else if (resname .eq. ' LI') then
               remark = 'HETATM'
            else if (resname .eq. '  F') then
               remark = 'HETATM'
            else if (resname .eq. ' NA') then
               remark = 'HETATM'
            else if (resname .eq. ' MG') then
               remark = 'HETATM'
            else if (resname .eq. ' CL') then
               remark = 'HETATM'
            else if (resname .eq. '  K') then
               remark = 'HETATM'
            else if (resname .eq. ' CA') then
               remark = 'HETATM'
            else if (resname .eq. ' FE') then
               remark = 'HETATM'
            else if (resname .eq. ' ZN') then
               remark = 'HETATM'
            else if (resname .eq. ' BR') then
               remark = 'HETATM'
            else if (resname .eq. '  I') then
               remark = 'HETATM'
            else if (residue.ne.reslast .or. resname.ne.namelast .or.
     &               chain.ne.chnlast .or. insert.ne.inslast) then
               nres = nres + 1
               reslast = residue
               namelast = resname
               chnlast = chain
               inslast = insert
               if (nres .gt. maxres) then
                  write (iout,70)  maxres
   70             format (/,' READPDB  --  The Maximum of',i6,
     &                       ' Residues has been Exceeded')
                  call fatal
               end if
               nseq = nres
               seq(nseq) = resname
            end if
            npdb = npdb + 1
            xpdb(npdb) = xx
            ypdb(npdb) = yy
            zpdb(npdb) = zz
            pdbrec(npdb) = remark
            pdbsym(npdb) = atmsymb
            pdbatm(npdb) = atmname
            pdbres(npdb) = resname
            pdbmod(npdb) = model
            resnum(npdb) = nres
            if (resname .eq. 'HOH')  resnum(npdb) = 0
            chnatm(npdb) = chain
   80       continue
         else if (remark .eq. 'HETATM') then
            next = 7
            call getnumb (record,serial,next)
            call gettext (record,atmsymb,next)
            call gettext (record,atmname,next)
            call gettext (record,altloc,next)
            call gettext (record,resname,next)
            call gettext (record,chain,next)
            call gettext (record,letter,next)
            call getnumb (record,residue,next)
            call gettext (record,insert,next)
            call getfloat (record,xx,next)
            call getfloat (record,yy,next)
            call getfloat (record,zz,next)
            call getfloat (record,occupy,next)
            call getfloat (record,bfac,next)
            call gettext (record,formal,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call getnumb (record,model,next)
            if (altloc .eq. '.')  altloc = ' '
            if (insert .eq. '?')  insert = ' '
            if (formal .eq. '?')  formal = ' '
            if (index(chnsym,chain) .eq. 0)  goto 90
            if (altloc.ne.' ' .and. altloc.ne.altsym)  goto 90
            if (insert.ne.' ' .and. index(instyp,insert).eq.0)  goto 90
            if (model .gt. 1)  goto 90
            call fixcif (resname,atmname)
            npdb = npdb + 1
            xpdb(npdb) = xx
            ypdb(npdb) = yy
            zpdb(npdb) = zz
            pdbrec(npdb) = remark
            pdbatm(npdb) = atmname
            pdbsym(npdb) = atmsymb
            pdbres(npdb) = resname
            pdbmod(npdb) = model
            resnum(npdb) = 0
            chnatm(npdb) = chain
   90       continue
         end if
         if (npdb .gt. maxatm) then
            write (iout,100)  maxatm
  100       format (/,' READCIF  --  The Maximum of',i6,
     &                 ' Atoms has been Exceeded')
            call fatal
         end if
      end do
  110 continue
c
c     set up crystal lattice if values were read from file
c
      if (xbx .gt. 1.0d0) then
         xbox = xbx
         ybox = ybx
         zbox = zbx
         alpha = aan
         beta = ban
         gamma = gan
         call unitcell
      end if
c
c     set the total sequence length and chain terminus sites
c
      if (npdb .ne. 0) then
         nchain = 0
         chnlast = '#'
         do i = 1, npdb
            if (pdbrec(i) .eq. 'ATOM  ') then
               letter = chnatm(i)
               if (letter .ne. chnlast) then
                  nchain = nchain + 1
                  ichain(1,nchain) = resnum(i)
                  chnnam(nchain) = letter
                  chnlast = letter
               else
                  ichain(2,nchain) = resnum(i)
               end if
            end if
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (chnatm)
c
c     find the type of species present in each chain
c
      do i = 1, nchain
         start = ichain(1,i)
         stop = ichain(2,i)
         chntyp(i) = 'GENERIC'
         do j = start, stop
            do k = 1, maxamino
               if (seq(j) .eq. amino(k)) then
                  chntyp(i) = 'PEPTIDE'
                  goto 120
               end if
            end do
            chntyp(i) = 'GENERIC'
            goto 130
  120       continue
         end do
  130    continue
         if (chntyp(i) .eq. 'GENERIC') then
            do j = start, stop
               do k = 1, maxnuc
                  if (seq(j) .eq. nuclz(k)) then
                     chntyp(i) = 'NUCLEIC'
                     goto 140
                  end if
               end do
               chntyp(i) = 'GENERIC'
               goto 150
  140          continue
            end do
  150       continue
         end if
      end do
c
c     get the three-letter sequence and code for each residue
c
      do i = 1, nchain
         start = ichain(1,i)
         stop = ichain(2,i)
         do j = start, stop
            do k = 1, maxamino
               if (seq(j) .eq. amino(k)) then
                  seqtyp(j) = k
                  goto 160
               end if
            end do
            do k = 1, maxnuc
               if (seq(j) .eq. nuclz(k)) then
                  seqtyp(j) = k
                  goto 160
               end if
            end do
            seq(j) = 'UNK'
            seqtyp(j) = 0
            if (chntyp(i) .eq. 'PEPTIDE')  seqtyp(j) = maxamino
            if (chntyp(i) .eq. 'NUCLEIC')  seqtyp(j) = maxnuc
  160       continue
         end do
      end do
c
c     set a pointer to the first and last atom of each residue
c
      nres = 0
      k = 0
      do i = 1, npdb
         if (pdbrec(i) .eq. 'ATOM  ') then
            if (resnum(i) .ne. k) then
               k = resnum(i)
               nres = nres + 1
               resatm(1,nres) = i
               if (nres .gt. 1)  resatm(2,nres-1) = i - 1
            end if
         end if
      end do
      if (nres .ge. 1)  resatm(2,nres) = npdb
c
c     close the input file and quit if no coordinates found
c
      if (.not. opened)  close (unit=icif)
      if (npdb .eq. 0)  abort = .true.
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine scancif  --  PDBx chains, alternates & inserts  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "scancif" reads the first model in a RCSB PDBx/mmCIF file
c     and finds chains, alternate sites, insertions and models
c
c
      subroutine scancif (icif)
      use iounit
      use pdb
      use sequen
      implicit none
      integer i,k,icif
      integer next,length
      integer nalt,nins
      integer model,modtemp
      logical exist,done
      character*1 letter
      character*1 chain,chnlast
      character*1 altloc,altlast
      character*1 insert,inslast
      character*6 remark
      character*20 blank,text
      character*20 chntemp
      character*20 alttyp
      character*20 instemp
      character*240 record
      character*240 string
c
c
c     initialize chain, alternate site and insertion lists
c
      nchain = 0
      nalt = 0
      nins = 0
      chnlast = '#'
      altlast = '#'
      inslast = '#'
      blank = '                    '
      chnsym = '####################'
      altsym = ' '
      alttyp = blank
      instyp = blank
      nmodel = 0
      imodel = 0
c
c     scan for multiple chains, alternate locations and inserts
c
      done = .false.
      do while (.not. done)
         read (icif,10,err=20,end=20)  record
   10    format (a240)
         remark = record(1:6)
         call upcase (remark)
         if (remark(1:5).eq.'ATOM ' .or. remark.eq.'HETATM') then
            next = 6
            if (remark .eq. 'HETATM')  next = 7
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,altloc,next)
            call gettext (record,letter,next)
            call gettext (record,chain,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,insert,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call gettext (record,letter,next)
            call getnumb (record,model,next)
            if (altloc .eq. '.')  altloc = ' '
            if (insert .eq. '?')  insert = ' '
            nmodel = max(model,nmodel)
            if (chain .ne. chnlast) then
               if (index(chnsym,chain) .eq. 0) then
                  nchain = nchain + 1
                  chnsym(nchain:nchain) = chain
                  chnlast = chain
               end if
            end if
            if (altloc .ne. altlast) then
               if (index(alttyp,altloc) .eq. 0) then
                  nalt = nalt + 1
                  alttyp(nalt:nalt) = altloc
                  altlast = altloc
               end if
            end if
            if (insert .ne. inslast) then
               if (index(instyp,insert) .eq. 0) then
                  nins = nins + 1
                  instyp(nins:nins) = insert
                  inslast = insert
               end if
            end if
         end if
      end do
   20 continue
      rewind (unit=icif)
c
c     find out which of the multiple chains will be used
c
      if (nchain .gt. 1) then
         call nextarg (chntemp,exist)
         if (.not. exist) then
            chntemp = blank
            if (chnsym(1:1) .eq. ' ') then
               string = 'BLANK=@'
               length = 7
            else
               string(1:1) = chnsym(1:1)
               length = 1
            end if
            do i = 2, nchain
               if (chnsym(i:i) .eq. ' ') then
                  string = string(1:length)//' BLANK=@'
                  length = length + 8
               else
                  string = string(1:length)//' '//chnsym(i:i)
                  length = length + 2
               end if
            end do
            string = string(1:length)//' [ALL]'
            length = length + 6
            write (iout,30)  string(1:length)
   30       format (/,' Enter the Chain Names to Include',
     &                 ' (',a,') :  ',$)
            read (input,40)  chntemp
   40       format (a20)
         end if
         call upcase (chntemp)
         next = 1
         call gettext (chntemp,text,next)
         if (text.eq.blank .or. text(1:3).eq.'ALL') then
            chnsym = chnsym(1:nchain)
         else
            do i = 1, nchain
               chain = chnsym(i:i)
               if (chain .eq. ' ')  chain = '@'
               k = index(chntemp,chain)
               if (k .eq. 0)  chnsym(i:i) = '#'
            end do
            chntemp = chnsym
            k = 0
            do i = 1, nchain
               chain = chntemp(i:i)
               if (chain .eq. '@')  chain = ' '
               if (chain .ne. '#') then
                  k = k + 1
                  chnsym(k:k) = chain
               end if
            end do
            nchain = k
         end if
      end if
      do i = nchain+1, 20
         chnsym(i:i) = '#'
      end do
c
c     find out which of alternate locations will be used
c
      if (nalt .gt. 0) then
         call nextarg (altsym,exist)
         if (.not. exist) then
            string(1:3) = '['//alttyp(1:1)//']'
            length = 3
            do i = 2, nalt
               string = string(1:length)//' '//alttyp(i:i)
               length = length + 2
            end do
            write (iout,50)  string(1:length)
   50       format (/,' Enter a Set of Alternate Atom Locations',
     &                 ' from (',a,') :  ',$)
            read (input,60)  record
   60       format (a240)
            next = 1
            call gettext (record,altsym,next)
         end if
         if (altsym .eq. ' ')  altsym = alttyp(1:1)
         call upcase (altsym)
      end if
c
c     find out which of the insert records will be used
c
      if (nins .gt. 0) then
         call nextarg (instemp,exist)
         if (.not. exist) then
            instemp = blank
            string(1:1) = instyp(1:1)
            length = 1
            do i = 2, nins
               string = string(1:length)//' '//instyp(i:i)
               length = length + 2
            end do
            string = string(1:length)//' [ALL] NONE'
            length = length + 11
            write (iout,70)  string(1:length)
   70       format (/,' Enter the Insert Records to Include',
     &                 ' (',a,') :  ',$)
            read (input,80)  instemp
   80       format (a20)
         end if
         call upcase (instemp)
         next = 1
         call gettext (instemp,text,next)
         if (text.eq.blank .or. text.eq.'ALL ') then
            instyp = instyp(1:nins)
         else if (text .eq. 'NONE ') then
            instyp = blank
         else
            instyp = instemp
         end if
      end if
c
c     find out which of the multiple models will be used
c
      if (nmodel .gt. 1) then
         call nextarg (string,exist)
         read (string,*,err=90,end=90)  modtemp
   90    continue
         if (.not. exist) then
            modtemp = 0
            write (iout,100)
  100       format (/,' Enter the Structural Model to Extract',
     &                 ' [0=All] :  ',$)
            read (input,110)  modtemp
  110       format (i10)
         end if
         if (modtemp .ne. 0)  nmodel = 1
         if (modtemp .eq. 0)  modtemp = 1
         imodel = modtemp
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fixcif  --  shift CIF atom names to PDB names  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fixcif" corrects CIF atom name entries by shifting them to
c     align with the standard PDB values
c
c
      subroutine fixcif (resname,atmname)
      implicit none
      character*1 first,last
      character*3 resname
      character*4 atmname
c
c
c     shift left-justified CIF names to standard PDB names
c
      first = atmname(1:1)
      if (first.ge.'A' .and.  first.le.'Z') then
         last = atmname(4:4)
         if (last .eq. ' ')  atmname = ' '//atmname(1:3)
      end if
c
c     convert unusual PDB names to their standard forms
c
      call fixpdb (resname,atmname)
      return
      end
