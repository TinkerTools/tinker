c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine readseq  --  read biopolymer sequence file  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "readseq" gets a biopolymer sequence containing one or more
c     separate chains from an external file; all lines containing
c     sequence must begin with the starting sequence number, the
c     actual sequence is read from subsequent nonblank characters
c
c
      subroutine readseq (iseq)
      implicit none
      include 'sizes.i'
      include 'files.i'
      include 'iounit.i'
      include 'resdue.i'
      include 'sequen.i'
      integer i,j,iseq
      integer length,number
      integer next,trimtext
      logical exist,opened,done
      character*1 letter
      character*3 word
      character*120 seqfile
      character*120 record
c
c
c     open the input file if it has not already been done
c
      inquire (unit=iseq,opened=opened)
      if (.not. opened) then
         seqfile = filename(1:leng)//'.seq'
         call version (seqfile,'old')
         inquire (file=seqfile,exist=exist)
         if (exist) then
            open (unit=iseq,file=seqfile,status='old')
            rewind (unit=iseq)
         else
            write (iout,10)
   10       format (/,' READSEQ  --  Unable to Find the Biopolymer',
     &                 ' Sequence File')
            call fatal
         end if
      end if
c
c     zero out the number and type of residues
c
      nseq = 0
      nchain = 0
      do i = 1, maxres
         seq(i) = '   '
      end do
c
c     read in the biopolymer sequence file
c
      do while (.true.)
         read (iseq,20,err=30,end=30)  record
   20    format (a120)
         length = trimtext (record)
         next = 1
         call gettext (record,letter,next)
         if (letter.ge.'0' .and. letter.le.'9') then
            next = 1
            letter = ' '
         end if
         call getnumb (record,number,next)
         if (number .eq. 1) then
            nchain = nchain + 1
            ichain(1,nchain) = nseq + 1
            chnnam(nchain) = letter
         end if
         done = .false.
         do while (.not. done)
            call getword (record,word,next)
            if (word .eq. '   ') then
               done = .true.
            else
               nseq = nseq + 1
               seq(nseq) = word
            end if
         end do
      end do
   30 continue
c
c     set the last residue in each sequence chain
c
      do i = 1, nchain-1
         ichain(2,i) = ichain(1,i+1) - 1
      end do
      if (nchain .ne. 0)  ichain(2,nchain) = nseq
c
c     find the residue type for each sequence element
c
      do i = 1, nseq
         seqtyp(i) = 0
         do j = 1, maxamino
            if (seq(i) .eq. amino(j)) then
               seqtyp(i) = j
               goto 40
            end if
         end do
         do j = 1, maxnuc
            if (seq(i) .eq. nuclz(j)) then
               seqtyp(i) = j
               goto 40
            end if
         end do
   40    continue
      end do
      if (.not. opened)  close (unit=iseq)
      return
      end
