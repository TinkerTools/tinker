c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  files.i  --  name and number of current structure files  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     nprior     number of previously existing cycle files
c     ldir       length in characters of the directory name
c     leng       length in characters of the base filename
c     filename   base filename used by default for all files
c     outfile    output filename used for intermediate results
c
c
      integer nprior,ldir,leng
      character*120 filename,outfile
      common /files/ nprior,ldir,leng,filename,outfile
