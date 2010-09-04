c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine readpdb  --  input of Protein Data Bank file  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readpdb" gets a set of Protein Data Bank coordinates
c     from an external disk file
c
c
      subroutine readpdb (ipdb)
      implicit none
      include 'sizes.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'pdb.i'
      include 'sequen.i'
      include 'titles.i'
      integer i,ipdb,nres
      integer index,serial
      integer trimtext
      integer residue,reslast
      real*8 xx,yy,zz
      logical exist,opened
      character*1 chain,chnlast
      character*1 altloc
      character*1 insert,inslast
      character*1 letter
      character*1 chnatm(maxatm)
      character*3 resname
      character*3 namelast
      character*4 atmname
      character*6 remark
      character*120 pdbfile
      character*120 record
      character*120 string
c
c
c     open the input file if it has not already been done
c
      inquire (unit=ipdb,opened=opened)
      if (.not. opened) then
         pdbfile = filename(1:leng)//'.pdb'
         call version (pdbfile,'old')
         inquire (file=pdbfile,exist=exist)
         if (exist) then
            open (unit=ipdb,file=pdbfile,status='old')
            rewind (unit=ipdb)
         else
            write (iout,10)
   10       format (/,' READPDB  --  Unable to Find the Protein',
     &                 ' Data Bank File')
            call fatal
         end if
      end if
c
c     get alternate sites, chains and insertions to be used
c
      call scanpdb (ipdb)
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
c     process individual atoms from the Protein Data Bank file
c
      do while (.true.)
         read (ipdb,20,err=90,end=90)  record
   20    format (a120)
         call upcase (record)
         remark = record(1:6)
         if (remark .eq. 'HEADER') then
            title = record(11:70)
            ltitle = trimtext (title)
         else if (remark .eq. 'TITLE ') then
            if (ltitle .eq. 0) then
               title = record(11:70)
               ltitle = trimtext (title)
            end if
         else if (remark .eq. 'ATOM  ') then
            string = record(7:120)
            read (string,30)  serial,atmname,altloc,resname,
     &                        chain,residue,insert,xx,yy,zz
   30       format (i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3)
            if (resname(1:2) .eq. '  ')  resname = resname(3:3)
            if (resname(1:1) .eq. ' ')  resname = resname(2:3)
            if (chain.ne.' ' .and. index(chntyp,chain).eq.0)  goto 50
            if (altloc.ne.' ' .and. altloc.ne.altsym)  goto 50
            if (insert.ne.' ' .and. index(instyp,insert).eq.0)  goto 50
            call fixpdb (resname,atmname)
            if (residue.ne.reslast .or. resname.ne.namelast .or.
     &            chain.ne.chnlast .or. insert.ne.inslast) then
               nres = nres + 1
               reslast = residue
               namelast = resname
               chnlast = chain
               inslast = insert
               if (nres .gt. maxres) then
                  write (iout,40)  maxres
   40             format (/,' READPDB  --  The Maximum of',i6,
     &                       ' Residues has been Exceeded')
                  call fatal
               end if
            end if
            npdb = npdb + 1
            xpdb(npdb) = xx
            ypdb(npdb) = yy
            zpdb(npdb) = zz
            pdbtyp(npdb) = remark
            atmnam(npdb) = atmname
            resnam(npdb) = resname
            resnum(npdb) = nres
            chnatm(npdb) = chain
   50       continue
         else if (remark .eq. 'HETATM') then
            string = record(7:120)
            read (string,60)  serial,atmname,altloc,resname,
     &                        chain,residue,insert,xx,yy,zz
   60       format (i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3)
            if (resname(1:2) .eq. '  ')  resname = resname(3:3)
            if (resname(1:1) .eq. ' ')  resname = resname(2:3)
            if (chain.ne.' ' .and. index(chntyp,chain).eq.0)  goto 70
            if (altloc.ne.' ' .and. altloc.ne.altsym)  goto 70
            if (insert.ne.' ' .and. index(instyp,insert).eq.0)  goto 70
            call fixpdb (resname,atmname)
            npdb = npdb + 1
            xpdb(npdb) = xx
            ypdb(npdb) = yy
            zpdb(npdb) = zz
            pdbtyp(npdb) = remark
            atmnam(npdb) = atmname
            resnam(npdb) = resname
            resnum(npdb) = 0
            chnatm(npdb) = chain
   70       continue
         else if (remark .eq. 'ENDMDL') then
            goto 90
         else if (remark .eq. 'END   ') then
            goto 90
         end if
         if (npdb .gt. maxatm) then
            write (iout,80)  maxatm
   80       format (/,' READPDB  --  The Maximum of',i6,
     &                 ' Atoms has been Exceeded')
            call fatal
         end if
      end do
   90 continue
c
c     set the total sequence length and chain termini information
c
      if (npdb .ne. 0) then
         nseq = npdb
         nchain = 0
         chnlast = '#'
         do i = 1, npdb
            if (pdbtyp(i) .eq. 'ATOM  ') then
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
c     close the PDB file and quit if there are no coordinates
c
      if (npdb .eq. 0)  abort = .true.
      if (.not. opened)  close (unit=ipdb)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine scanpdb  --  PDB chains, alternates and inserts  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "scanpdb" reads the first model in a Protein Data Bank file and
c     sets chains, alternate sites and insertion records to be used
c
c
      subroutine scanpdb (ipdb)
      implicit none
      include 'sizes.i'
      include 'iounit.i'
      include 'pdb.i'
      include 'sequen.i'
      integer i,ipdb
      integer length,next
      integer nalt,nins
      logical exist,done,first
      character*1 chain,chnlast
      character*1 altloc,altlast
      character*1 insert,inslast
      character*6 remark
      character*20 blank,text
      character*20 chntemp
      character*20 alttyp
      character*20 instemp
      character*120 record
      character*120 string
      save first
      data first  / .true. /
c
c
c     only proceed if this routine was not already used
c
      if (.not. first)  return
      first = .false.
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
      chntyp = '####################'
      altsym = ' '
      alttyp = blank
      instyp = blank
c
c     scan for multiple chains, alternate locations and inserts
c
      done = .false.
      do while (.not. done)
         read (ipdb,10,err=30,end=30)  record
   10    format (a120)
         call upcase (record)
         remark = record(1:6)
         string = record(7:120)
         if (remark.eq.'ATOM  ' .or. remark.eq.'HETATM') then
            read (string,20)  altloc,chain,insert
   20       format (10x,a1,4x,a1,4x,a1)
            if (chain .ne. chnlast) then
               if (index(chntyp,chain) .eq. 0) then
                  nchain = nchain + 1
                  chntyp(nchain:nchain) = chain
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
         else if (remark .eq. 'ENDMDL') then
            done = .true.
         else if (remark .eq. 'END   ') then
            done = .true.
         end if
      end do
   30 continue
      rewind (unit=ipdb)
c
c     find out which of the multiple chains will be used
c
      if (nchain .gt. 1) then
         call nextarg (chntemp,exist)
         if (.not. exist) then
            chntemp = blank
            if (chntyp(1:1) .eq. ' ') then
               string = 'BLANK'
               length = 5
            else
               string(1:1) = chntyp(1:1)
               length = 1
            end if
            do i = 2, nchain
               if (chntyp(i:i) .eq. ' ') then
                  string = string(1:length)//' BLANK'
                  length = length + 6
               else
                  string = string(1:length)//' '//chntyp(i:i)
                  length = length + 2
               end if
            end do
            string = string(1:length)//' [ALL]'
            length = length + 6
            write (iout,40)  string(1:length)
   40       format (/,' Enter the Chain Names to Include',
     &                 ' (',a,') :  ',$)
            read (input,50)  chntemp
   50       format (a20)
         end if
         call upcase (chntemp)
         next = 1
         call gettext (chntemp,text,next)
         if (text.eq.blank .or. text.eq.'ALL ') then
            chntyp = chntyp(1:nchain)
         else
            nchain = 1
            chntyp = chntemp(1:1)
         end if
      end if
c
c     find out which of the alternate locations will be used
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
            write (iout,60)  string(1:length)
   60       format (/,' Enter a Set of Alternate Atom Locations',
     &                 ' from (',a,') :  ',$)
            read (input,70)  record
   70       format (a120)
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
            write (iout,80)  string(1:length)
   80       format (/,' Enter the Insert Records to Include',
     &                 ' (',a,') :  ',$)
            read (input,90)  instemp
   90       format (a20)
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
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine fixpdb  --  standard PDB atom and residue names  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "fixpdb" corrects problems with PDB files by converting residue
c     and atom names to the standard forms used by TINKER
c
c
      subroutine fixpdb (resname,atmname)
      implicit none
      include 'sizes.i'
      include 'resdue.i'
      integer i
      character*3 resname
      character*4 atmname
      character*7 restype
c
c
c     convert traditional 3-letter base names to PDB names
c
      if (resname .eq. 'ADE')  resname = 'A  '
      if (resname .eq. 'GUA')  resname = 'G  '
      if (resname .eq. 'CYT')  resname = 'C  '
      if (resname .eq. 'URA')  resname = 'U  '
      if (resname .eq. 'DAD')  resname = 'DA '
      if (resname .eq. 'DGU')  resname = 'DG '
      if (resname .eq. 'DCY')  resname = 'DC '
      if (resname .eq. 'THY')  resname = 'DT '
c
c     convert unusual names for protonated histidine residues
c
      if (resname .eq. 'HSD')  resname = 'HID'
      if (resname .eq. 'HSE')  resname = 'HIE'
      if (resname .eq. 'HSP')  resname = 'HIS'
      if (resname .eq. 'HSH')  resname = 'HIS'
      if (resname .eq. 'HIP')  resname = 'HIS'
      if (resname .eq. 'HIH')  resname = 'HIS'
c
c     convert unusual names for terminal capping residues
c
      if (resname .eq. 'NMA')  resname = 'NME'
c
c     convert nonstandard names for water molecules
c
      if (resname .eq. 'H2O')  resname = 'HOH'
      if (resname .eq. 'WAT')  resname = 'HOH'
c
c     decide whether residue is protein or nucleic acid
c
      restype = 'UNKNOWN'
      do i = 1, maxamino
         if (resname .eq. amino(i))  restype = 'PROTEIN'
      end do
      do i = 1, maxnuc
         if (resname .eq. nuclz(i))  restype = 'NUCLEIC'
      end do
c
c     convert any generically used unusual atom names
c
      if (atmname .eq. ' HN ')  atmname = ' H  '
c
c     convert unusual names in protein terminal residues
c
      if (restype .eq. 'PROTEIN') then
         if (atmname .eq. '1H  ')  atmname = ' H1 '
         if (atmname .eq. ' HN1')  atmname = ' H1 '
         if (atmname .eq. ' HT1')  atmname = ' H1 '
         if (atmname .eq. '2H  ')  atmname = ' H2 '
         if (atmname .eq. ' HN2')  atmname = ' H2 '
         if (atmname .eq. ' HT2')  atmname = ' H2 '
         if (atmname .eq. '3H  ')  atmname = ' H3 '
         if (atmname .eq. ' HN3')  atmname = ' H3 '
         if (atmname .eq. ' HT3')  atmname = ' H3 '
         if (atmname .eq. ' O1 ')  atmname = ' O  '
         if (atmname .eq. ' OT1')  atmname = ' O  '
         if (atmname .eq. 'OCT1')  atmname = ' O  '
         if (atmname .eq. ' O2 ')  atmname = ' OXT'
         if (atmname .eq. ' OT2')  atmname = ' OXT'
         if (atmname .eq. 'OCT2')  atmname = ' OXT'
         if (atmname .eq. ' OT ')  atmname = ' OXT'
      end if
c
c     convert unusual names common to many nucleotides
c
      if (restype .eq. 'NUCLEIC') then
         if (atmname .eq. ' O1P')  atmname = ' OP1'
         if (atmname .eq. ' O2P')  atmname = ' OP2'
         if (atmname .eq. ' O3P')  atmname = ' OP3'
         if (atmname .eq. '2HOP')  atmname = 'HOP2'
         if (atmname .eq. '3HOP')  atmname = 'HOP3'
         if (atmname .eq. ' C1*')  atmname = ' C1'''
         if (atmname .eq. ' C2*')  atmname = ' C2'''
         if (atmname .eq. ' C3*')  atmname = ' C3'''
         if (atmname .eq. ' C4*')  atmname = ' C4'''
         if (atmname .eq. ' C5*')  atmname = ' C5'''
         if (atmname .eq. ' O2*')  atmname = ' O2'''
         if (atmname .eq. ' O3*')  atmname = ' O3'''
         if (atmname .eq. ' O4*')  atmname = ' O4'''
         if (atmname .eq. ' O5*')  atmname = ' O5'''
         if (atmname .eq. ' H1*')  atmname = ' H1'''
         if (atmname .eq. ' H2*')  atmname = ' H2'''
         if (atmname .eq. '1H2*')  atmname = ' H2'''
         if (atmname .eq. '2H2*')  atmname = 'H2'''''
         if (atmname .eq. ' H3*')  atmname = ' H3'''
         if (atmname .eq. ' H4*')  atmname = ' H4'''
         if (atmname .eq. '1H5*')  atmname = ' H5'''
         if (atmname .eq. '2H5*')  atmname = 'H5'''''
         if (atmname .eq. '2HO*')  atmname = 'HO2'''
         if (atmname .eq. ' H3T')  atmname = 'HO3'''
      end if
c
c     glycine residue  (GLY)
c
      if (resname .eq. 'GLY') then
         if (atmname .eq. '1HA ')  atmname = ' HA2'
         if (atmname .eq. ' HA1')  atmname = ' HA3'
         if (atmname .eq. '2HA ')  atmname = ' HA3'
c
c     alanine residue  (ALA)
c
      else if (resname .eq. 'ALA') then
         if (atmname .eq. '1HB ')  atmname = ' HB1'
         if (atmname .eq. '2HB ')  atmname = ' HB2'
         if (atmname .eq. '3HB ')  atmname = ' HB3'
c
c     valine residue  (VAL)
c
      else if (resname .eq. 'VAL') then
         if (atmname .eq. '1HG1')  atmname = 'HG11'
         if (atmname .eq. '2HG1')  atmname = 'HG12'
         if (atmname .eq. '3HG1')  atmname = 'HG13'
         if (atmname .eq. '1HG2')  atmname = 'HG21'
         if (atmname .eq. '2HG2')  atmname = 'HG22'
         if (atmname .eq. '3HG2')  atmname = 'HG23'
c
c     leucine residue  (LEU)
c
      else if (resname .eq. 'LEU') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. '1HD1')  atmname = 'HD11'
         if (atmname .eq. '2HD1')  atmname = 'HD12'
         if (atmname .eq. '3HD1')  atmname = 'HD13'
         if (atmname .eq. '1HD2')  atmname = 'HD21'
         if (atmname .eq. '2HD2')  atmname = 'HD22'
         if (atmname .eq. '3HD2')  atmname = 'HD23'
c
c     isoleucine residue  (ILE)
c
      else if (resname .eq. 'ILE') then
         if (atmname .eq. ' CD ')  atmname = ' CD1'
         if (atmname .eq. '1HG1')  atmname = 'HG12'
         if (atmname .eq. 'HG11')  atmname = 'HG13'
         if (atmname .eq. '2HG1')  atmname = 'HG13'
         if (atmname .eq. '1HG2')  atmname = 'HG21'
         if (atmname .eq. '2HG2')  atmname = 'HG22'
         if (atmname .eq. '3HG2')  atmname = 'HG23'
         if (atmname .eq. '1HD1')  atmname = 'HD11'
         if (atmname .eq. ' HD1')  atmname = 'HD11'
         if (atmname .eq. '2HD1')  atmname = 'HD12'
         if (atmname .eq. ' HD2')  atmname = 'HD12'
         if (atmname .eq. '3HD1')  atmname = 'HD13'
         if (atmname .eq. ' HD3')  atmname = 'HD13'
c
c     serine residue  (SER)
c
      else if (resname .eq. 'SER') then
         if (atmname .eq. ' OG1')  atmname = ' OG '
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. ' HG1')  atmname = ' HG '
         if (atmname .eq. ' HOG')  atmname = ' HG '
c
c     threonine residue  (THR)
c
      else if (resname .eq. 'THR') then
         if (atmname .eq. ' OG ')  atmname = ' OG1'
         if (atmname .eq. ' CG ')  atmname = ' CG2'
         if (atmname .eq. ' HOG')  atmname = ' HG1'
         if (atmname .eq. 'HOG1')  atmname = ' HG1'
         if (atmname .eq. '1HG2')  atmname = 'HG21'
         if (atmname .eq. '2HG2')  atmname = 'HG22'
         if (atmname .eq. '3HG2')  atmname = 'HG23'
c
c     cysteine residue  (CYS)
c
      else if (resname .eq. 'CYS') then
         if (atmname .eq. ' SG1')  atmname = ' SG '
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. ' HG1')  atmname = ' HG '
         if (atmname .eq. ' HSG')  atmname = ' HG '
c
c     proline residue  (PRO)
c
      else if (resname .eq. 'PRO') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. '1HG ')  atmname = ' HG2'
         if (atmname .eq. ' HG1')  atmname = ' HG3'
         if (atmname .eq. '2HG ')  atmname = ' HG3'
         if (atmname .eq. '1HD ')  atmname = ' HD2'
         if (atmname .eq. ' HD1')  atmname = ' HD3'
         if (atmname .eq. '2HD ')  atmname = ' HD3'
         if (atmname .eq. ' HT1')  atmname = ' H2 '
         if (atmname .eq. ' HT2')  atmname = ' H3 '
c
c     phenylalanine residue  (PHE)
c
      else if (resname .eq. 'PHE') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
c
c     tyrosine residue  (TYR)
c
      else if (resname .eq. 'TYR') then
         if (atmname .eq. ' HOH')  atmname = ' HH '
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
c
c     tryptophan residue  (TRP)
c
      else if (resname .eq. 'TRP') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. ' HNE')  atmname = ' HE1'
c
c     histidine (HD and HE) residue  (HIS)
c
      else if (resname .eq. 'HIS') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. ' HD ')  atmname = ' HD2'
         if (atmname .eq. ' HE ')  atmname = ' HE1'
         if (atmname .eq. ' HND')  atmname = ' HD1'
         if (atmname .eq. 'HND1')  atmname = ' HD1'
         if (atmname .eq. ' HNE')  atmname = ' HE2'
         if (atmname .eq. 'HNE2')  atmname = ' HE2'
c
c     histidine (HD only) residue  (HID)
c
      else if (resname .eq. 'HID') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. ' HD ')  atmname = ' HD2'
         if (atmname .eq. ' HE ')  atmname = ' HE1'
         if (atmname .eq. ' HND')  atmname = ' HD1'
         if (atmname .eq. 'HND1')  atmname = ' HD1'
c
c     histidine (HE only) residue  (HIE)
c
      else if (resname .eq. 'HIE') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. ' HD ')  atmname = ' HD2'
         if (atmname .eq. ' HE ')  atmname = ' HE1'
         if (atmname .eq. ' HNE')  atmname = ' HE2'
         if (atmname .eq. 'HNE2')  atmname = ' HE2'
c
c     aspartic acid residue  (ASP)
c
      else if (resname .eq. 'ASP') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
c
c     asparagine residue  (ASN)
c
      else if (resname .eq. 'ASN') then
         if (atmname .eq. ' OD ')  atmname = ' OD1'
         if (atmname .eq. ' ND ')  atmname = ' ND2'
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. '1HD2')  atmname = 'HD21'
         if (atmname .eq. 'HND1')  atmname = 'HD21'
         if (atmname .eq. '2HD2')  atmname = 'HD22'
         if (atmname .eq. 'HND2')  atmname = 'HD22'
c
c     glutamic acid residue  (GLU)
c
      else if (resname .eq. 'GLU') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. '1HG ')  atmname = ' HG2'
         if (atmname .eq. ' HG1')  atmname = ' HG3'
         if (atmname .eq. '2HG ')  atmname = ' HG3'
c
c     glutamine residue  (GLN)
c
      else if (resname .eq. 'GLN') then
         if (atmname .eq. ' OE ')  atmname = ' OE1'
         if (atmname .eq. ' NE ')  atmname = ' NE2'
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. '1HG ')  atmname = ' HG2'
         if (atmname .eq. ' HG1')  atmname = ' HG3'
         if (atmname .eq. '2HG ')  atmname = ' HG3'
         if (atmname .eq. '1HE2')  atmname = 'HE21'
         if (atmname .eq. 'HNE1')  atmname = 'HE21'
         if (atmname .eq. '2HE2')  atmname = 'HE22'
         if (atmname .eq. 'HNE2')  atmname = 'HE22'
c
c     methionine residue  (MET)
c
      else if (resname .eq. 'MET') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. '1HG ')  atmname = ' HG2'
         if (atmname .eq. ' HG1')  atmname = ' HG3'
         if (atmname .eq. '2HG ')  atmname = ' HG3'
         if (atmname .eq. '1HE ')  atmname = ' HE1'
         if (atmname .eq. '2HE ')  atmname = ' HE2'
         if (atmname .eq. '3HE ')  atmname = ' HE3'
c
c     lysine residue  (LYS)
c
      else if (resname .eq. 'LYS') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. '1HG ')  atmname = ' HG2'
         if (atmname .eq. ' HG1')  atmname = ' HG3'
         if (atmname .eq. '2HG ')  atmname = ' HG3'
         if (atmname .eq. '1HD ')  atmname = ' HD2'
         if (atmname .eq. ' HD1')  atmname = ' HD3'
         if (atmname .eq. '2HD ')  atmname = ' HD3'
         if (atmname .eq. '1HE ')  atmname = ' HE2'
         if (atmname .eq. ' HE1')  atmname = ' HE3'
         if (atmname .eq. '2HE ')  atmname = ' HE3'
         if (atmname .eq. '1HZ ')  atmname = ' HZ1'
         if (atmname .eq. 'HNZ1')  atmname = ' HZ1'
         if (atmname .eq. '2HZ ')  atmname = ' HZ2'
         if (atmname .eq. 'HNZ2')  atmname = ' HZ2'
         if (atmname .eq. '3HZ ')  atmname = ' HZ3'
         if (atmname .eq. 'HNZ3')  atmname = ' HZ3'
c
c     arginine residue  (ARG)
c
      else if (resname .eq. 'ARG') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. '1HG ')  atmname = ' HG2'
         if (atmname .eq. ' HG1')  atmname = ' HG3'
         if (atmname .eq. '2HG ')  atmname = ' HG3'
         if (atmname .eq. '1HD ')  atmname = ' HD2'
         if (atmname .eq. ' HD1')  atmname = ' HD3'
         if (atmname .eq. '2HD ')  atmname = ' HD3'
         if (atmname .eq. '1HH1')  atmname = 'HH11'
         if (atmname .eq. 'HN11')  atmname = 'HH11'
         if (atmname .eq. '2HH1')  atmname = 'HH12'
         if (atmname .eq. 'HN12')  atmname = 'HH12'
         if (atmname .eq. '1HH2')  atmname = 'HH21'
         if (atmname .eq. 'HN21')  atmname = 'HH21'
         if (atmname .eq. '2HH2')  atmname = 'HH22'
         if (atmname .eq. 'HN22')  atmname = 'HH22'
c
c     ornithine residue  (ORN)
c
      else if (resname .eq. 'ORN') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. '1HG ')  atmname = ' HG2'
         if (atmname .eq. ' HG1')  atmname = ' HG3'
         if (atmname .eq. '2HG ')  atmname = ' HG3'
         if (atmname .eq. '1HD ')  atmname = ' HD2'
         if (atmname .eq. ' HD1')  atmname = ' HD3'
         if (atmname .eq. '2HD ')  atmname = ' HD3'
         if (atmname .eq. '1HE ')  atmname = ' HE1'
         if (atmname .eq. 'HNE1')  atmname = ' HE1'
         if (atmname .eq. '2HE ')  atmname = ' HE2'
         if (atmname .eq. 'HNE2')  atmname = ' HE2'
         if (atmname .eq. '3HE ')  atmname = ' HE3'
         if (atmname .eq. 'HNE3')  atmname = ' HE3'
c
c     methylalanine residue  (AIB)
c
      else if (resname .eq. 'AIB') then
         if (atmname .eq. '1HB1')  atmname = 'HB11'
         if (atmname .eq. '2HB1')  atmname = 'HB12'
         if (atmname .eq. '3HB1')  atmname = 'HB13'
         if (atmname .eq. '1HB2')  atmname = 'HB21'
         if (atmname .eq. '2HB2')  atmname = 'HB22'
         if (atmname .eq. '3HB2')  atmname = 'HB23'
c
c     pyroglutamic acid residue  (PCA)
c
      else if (resname .eq. 'PCA') then
         if (atmname .eq. '1HB ')  atmname = ' HB2'
         if (atmname .eq. ' HB1')  atmname = ' HB3'
         if (atmname .eq. '2HB ')  atmname = ' HB3'
         if (atmname .eq. '1HG ')  atmname = ' HG2'
         if (atmname .eq. ' HG1')  atmname = ' HG3'
         if (atmname .eq. '2HG ')  atmname = ' HG3'
c
c     N-terminal acetyl residue  (ACE)
c
      else if (resname .eq. 'ACE') then
         if (atmname .eq. ' CY ')  atmname = ' C  '
         if (atmname .eq. ' CAY')  atmname = ' CH3'
         if (atmname .eq. ' CA ')  atmname = ' CH3'
         if (atmname .eq. ' OY ')  atmname = ' O  '
         if (atmname .eq. '1H  ')  atmname = ' H1 '
         if (atmname .eq. ' HY1')  atmname = ' H1 '
         if (atmname .eq. 'HH31')  atmname = ' H1 '
         if (atmname .eq. '2H  ')  atmname = ' H2 '
         if (atmname .eq. ' HY2')  atmname = ' H2 '
         if (atmname .eq. 'HH32')  atmname = ' H2 '
         if (atmname .eq. '3H  ')  atmname = ' H3 '
         if (atmname .eq. ' HY3')  atmname = ' H3 '
         if (atmname .eq. 'HH33')  atmname = ' H3 '
c
c     N-terminal formyl residue  (FOR)
c
      else if (resname .eq. 'FOR') then
         if (atmname .eq. ' CY ')  atmname = ' C  '
         if (atmname .eq. ' OY ')  atmname = ' O  '
         if (atmname .eq. ' HY ')  atmname = ' H  '
c
c     C-terminal N-methylamide residue  (NME)
c
      else if (resname .eq. 'NME') then
         if (atmname .eq. ' NT ')  atmname = ' N  '
         if (atmname .eq. ' CT ')  atmname = ' CH3'
         if (atmname .eq. ' CAT')  atmname = ' CH3'
         if (atmname .eq. ' CA ')  atmname = ' CH3'
         if (atmname .eq. ' HNT')  atmname = ' H  '
         if (atmname .eq. '1H  ')  atmname = ' H1 '
         if (atmname .eq. '1HA ')  atmname = ' H1 '
         if (atmname .eq. ' HT1')  atmname = ' H1 '
         if (atmname .eq. 'HH31')  atmname = ' H1 '
         if (atmname .eq. '2H  ')  atmname = ' H2 '
         if (atmname .eq. '2HA ')  atmname = ' H2 '
         if (atmname .eq. ' HT2')  atmname = ' H2 '
         if (atmname .eq. 'HH32')  atmname = ' H2 '
         if (atmname .eq. '3H  ')  atmname = ' H3 '
         if (atmname .eq. '3HA ')  atmname = ' H3 '
         if (atmname .eq. ' HT3')  atmname = ' H3 '
         if (atmname .eq. 'HH33')  atmname = ' H3 '
c
c     C-terminal amide residue  (NH2)
c
      else if (resname .eq. 'NH2') then
         if (atmname .eq. ' NT ')  atmname = ' N  '
         if (atmname .eq. '1H  ')  atmname = ' H1 '
         if (atmname .eq. '2H  ')  atmname = ' H2 '
         if (atmname .eq. ' HT1')  atmname = ' H1 '
         if (atmname .eq. ' HT2')  atmname = ' H2 '
c
c     adenosine residue  (A)
c
      else if (resname .eq. 'A  ') then
         if (atmname .eq. '1H6 ')  atmname = ' H61'
         if (atmname .eq. '2H6 ')  atmname = ' H62'
c
c     guanosine residue  (G)
c
      else if (resname .eq. 'G  ') then
         if (atmname .eq. '1H2 ')  atmname = ' H21'
         if (atmname .eq. '2H2 ')  atmname = ' H22'
c
c     cytidine residue  (C)
c
      else if (resname .eq. 'C  ') then
         if (atmname .eq. '1H4 ')  atmname = ' H41'
         if (atmname .eq. '2H4 ')  atmname = ' H42'
c
c     deoxyadenosine residue  (DA)
c
      else if (resname .eq. 'DA ') then
         if (atmname .eq. '1H6 ')  atmname = ' H61'
         if (atmname .eq. '2H6 ')  atmname = ' H62'
c
c     deoxyguanosine residue  (DG)
c
      else if (resname .eq. 'DG ') then
         if (atmname .eq. '1H2 ')  atmname = ' H21'
         if (atmname .eq. '2H2 ')  atmname = ' H22'
c
c     deoxycytidine residue  (DC)
c
      else if (resname .eq. 'DC ') then
         if (atmname .eq. '1H4 ')  atmname = ' H41'
         if (atmname .eq. '2H4 ')  atmname = ' H42'
c
c     deoxythymidine residue  (DT)
c
      else if (resname .eq. 'DT ') then
         if (atmname .eq. ' C5M')  atmname = ' C7 '
         if (atmname .eq. '1H5M')  atmname = ' H71'
         if (atmname .eq. '2H5M')  atmname = ' H72'
         if (atmname .eq. '3H5M')  atmname = ' H73'
      end if
      return
      end
