c
c
c     #############################################################
c     ##                  COPYRIGHT (C) 1999 by                  ##
c     ##  Marina A. Vorobieva, Nina N. Sokolova & Jay W. Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  nucleo.i  --  parameters for nucleic acid structure  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     bkbone    phosphate backbone angles for each nucleotide
c     glyco     glycosidic torsional angle for each nucleotide
c     pucker    sugar pucker, either 2=2'-endo or 3=3'-endo
c     dblhlx    flag to mark system as nucleic acid double helix
c     deoxy     flag to mark deoxyribose or ribose sugar units
c     hlxform   helix form (A, B or Z) of polynucleotide strands
c
c
      integer pucker
      real*8 bkbone,glyco
      logical dblhlx,deoxy
      character*1 hlxform
      common /nucleo/ bkbone(6,maxres),glyco(maxres),pucker(maxres),
     &                dblhlx,deoxy(maxres),hlxform
