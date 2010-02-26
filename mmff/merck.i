c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2007 by Nicolas Staelens & Jay W. Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  merck.i  --  parameters specific for MMFF force field  ##
c     ##                                                         ##
c     #############################################################
c
c
c     BT_1     atom pairs having a Bond Type = 1
c     nlignes  number of atom pairs having a Bond Type = 1
c     eqclass  table containing the atom class equivalencies
c              to find 'default' parameters when angle bending,
c              out-of-plane, or torsion parameters are missing
c              (see T. A. Halgren, J. Comput. Chem., 17,
c               490-519, '95, Table IV, page 514)
c
c     crd      number of attached neighbours   |
c     val      "valency"                       |  see T. A. Halgren,
c     pilp     if 0 : no lone pair             |  J. Comput. Chem,
c              if 1 : one or more lone pair(s) |  17, 616-645 (1995)
c     mltb     multibond indicator             |
c     arom     aromaticity indicator           |
c     lin      linearity indicator             |
c     sbmb     single-bond/multiple-bond       |
c              indicator                       |
c
c     MMFFAROM     aromatic rings parameters
c     MMFFAROMCAT  cationic aromatic rings parameters
c     MMFFAROMAN   anionic aromatic rings parameters
c
c     MMFF_kb   idem, as a table (the force constant is at the 
c                  intersection of the atom classes of the two atoms)
c     MMFF_kb1  idem, with Bond Type = 1
c     MMFF_b0   idem, as a table (the bond length parameter is at the
c                  intersection of the atom classes of the two atoms)
c     MMFF_b01  idem, with Bond Type = 1
c
c     rad0    covalent radius of atom             |
c             (R. Blom and A. Haaland, J. Mol.    |  used with MMFF
c             Struct., 128, 21-27 '85)            |  empirical rules
c     paulel  Pauling-scaled electronegativities  |  for building
c             as defined by Allred and Rochow     |  missing bond
c     r0ref   reference bond length               |  stretch values
c     kbref   reference force constant            |
c
c     MMFF_ka   idem with an entry for the class of each of the
c                atoms of the angle
c     MMFF_ka1  idem with one of the two bonds having a Bond Type = 1
c     MMFF_ka2  idem with the two bonds having Bond Types of 1;
c                 the sum is 2
c     MMFF_ka3  idem with both bonds having a Bond Type = 0
c                 and the atoms belonging to a 3-membered ring
c     MMFF_ka4  idem with both bonds having a Bond Type = 0
c                 and the atoms belonging to a 4-membered ring
c     MMFF_ka5  idem with one of the bonds having a Bond Type = 1
c                 and the atoms belonging to a 3-membered ring
c     MMFF_ka6  idem with both bonds having a Bond Type = 1
c                 and the atoms belonging to a 3-membered ring
c     MMFF_ka7  idem with one of the bonds having a Bond Type = 1
c                 and the atoms belonging to a 4-membered ring
c     MMFF_ka8  idem with both bonds having a Bond Type = 1
c                 and the atoms belonging to a 4-membered ring
c     MMFF_teta0   idem with an entry for the class of each of the
c                    atoms of the angle
c     MMFF_teta01  idem with one of the bonds having a Bond Type = 1
c     MMFF_teta02  idem with both bonds having Bond Types of 1;
c                    the sum is 2
c     MMFF_teta03  idem with both bonds having Bond Types of 0
c                    and forming a 3-membered ring
c     MMFF_teta04  idem with both bonds having Bond Types of 0
c                    and forming a 4-membered ring
c     MMFF_teta05  idem with one of the bonds having Bond Types of 1
c                    and forming a 3-membered ring
c     MMFF_teta06  idem with both bonds having Bond Types of 1
c                    and forming a 3-membered ring
c     MMFF_teta07  idem with one of the bonds having Bond Types of 1
c                    and forming a 4-membered ring
c     MMFF_teta08  idem with both bonds having Bond Types of 1
c
c     stbn_abc   stretch-bending parameters for each atom class (A-B-C)
c     stbn_cba   stretch-bending parameters for each atom class (C-B-A)
c     with A-B having a Bond Type = 1, Stretch-Bend Type = 1,
c     stbn_abc1   stretch-bending parameters for each atom class (A-B-C)
c     stbn_cba1   stretch-bending parameters for each atom class (C-B-A)
c     with B-C having a Bond Type = 1, Stretch-Bend Type = 2
c     stbn_abc2   stretch-bending parameters for each atom class (A-B-C)
c     stbn_cba2   stretch-bending parameters for each atom class (C-B-A)
c     with A-B and B-C having a Bond Type = 1, Stretch-Bend Type = 3
c     stbn_abc3   stretch-bending parameters for each atom class (A-B-C)
c     stbn_cba3   stretch-bending parameters for each atom class (C-B-A)
c     with both Bond Type = 0 and 3-membered ring, 
c     Stretch-Bend Type = 5
c     stbn_abc5   stretch-bending parameters for each atom class (A-B-C)
c     stbn_cba5   stretch-bending parameters for each atom class (C-B-A)
c     with both Bond Type = 0 and 4-membered ring,         
c     Stretch-Bend Type = 4
c     stbn_abc4   stretch-bending parameters for each atom class (A-B-C)
c     stbn_cba4   stretch-bending parameters for each atom class (C-B-A)
c     with A-B having a Bond Type = 1 and 3-membered ring,
c     Stretch-Bend Type = 6
c     stbn_abc6   stretch-bending parameters for each atom class (A-B-C)
c     stbn_cba6   stretch-bending parameters for each atom class (C-B-A)
c     with B-C having a Bond Type = 1 and 3-membered ring,
c     Stretch-Bend Type = 7
c     stbn_abc7   stretch-bending parameters for each atom class (A-B-C)
c     stbn_cba7   stretch-bending parameters for each atom class (C-B-A)
c     with both Bond Type = 1 and 3-membered ring,         
c     Stretch-Bend Type = 8
c     stbn_abc8   stretch-bending parameters for each atom class (A-B-C)
c     stbn_cba8   stretch-bending parameters for each atom class (C-B-A)
c     with A-B having a Bond Type = 1 and 4-membered ring,
c     Stretch-Bend Type = 9
c     stbn_abc9   stretch-bending parameters for each atom class (A-B-C)
c     stbn_cba9   stretch-bending parameters for each atom class (C-B-A)
c     with B-C having a Bond Type = 1 and 4-membered ring,
c     Stretch-Bend Type = 10
c     stbn_abc10   stretch-bending parameters for each atom class (A-B-C)
c     stbn_cba10   stretch-bending parameters for each atom class (C-B-A)
c     with both Bond Type = 1 and 4-membered ring,
c     Stretch-Bend Type = 11
c     stbn_abc11   stretch-bending parameters for each atom class (A-B-C)
c     stbn_cba11   stretch-bending parameters for each atom class (C-B-A)
c
c     t1_1     torsional parameters for 1-fold, MMFF Torsion Type = 1
c     t1_2     torsional parameters for 1-fold, MMFF Torsion Type = 2
c     t2_1     torsional parameters for 2-fold, MMFF Torsion Type = 1
c     t2_2     torsional parameters for 2-fold, MMFF Torsion Type = 2
c     t3_1     torsional parameters for 3-fold, MMFF Torsion Type = 1
c     t3_2     torsional parameters for 3-fold, MMFF Torsion Type = 2
c     kt_1     string of classes for torsions, MMFF Torsion Type = 1
c     kt_2     string of classes for torsions, MMFF Torsion Type = 2
c
c     G        scale factors for calculation of MMFF eps
c     alph     atomic polarizabilities for calculation of MMFF eps
c     Nn       effective numbers of valence electrons 
c                for calculation of MMFF eps 
c     DA       donor/acceptor atom classes
c
c     bci      bond charge increments for the building of atom charges
c                 in MMFF depending on the attached neighbours 
c     bci_1    bond charge increments when Bond Type = 1   
c     pbci     partial BCI for building missing BCI's
c     fcadj    formal charge adjustment factor
c
c
      integer BT_1
      integer nlignes
      integer eqclass
      integer crd,val,pilp,mltb
      integer arom,lin,sbmb
      integer MMFFAROM,MMFFAROMCAT
      integer MMFFAROMAN
      common /merck1/ BT_1(500,2),nlignes,eqclass(500,5),crd(100),
     &                val(100),pilp(100),mltb(100),arom(100),lin(100),
     &                sbmb(100),MMFFAROM(maxtyp,6),
     &                MMFFAROMCAT(maxtyp,6),MMFFAROMAN(maxtyp,6)
c
      real*8 MMFF_kb,MMFF_kb1
      real*8 MMFF_b0,MMFF_b01
      real*8 rad0,r0ref
      real*8 kbref,paulel
      common /merck2/ MMFF_kb(100,100),MMFF_kb1(100,100),
     &                MMFF_b0(100,100),MMFF_b01(100,100),
     &                rad0(120),r0ref(53,53),kbref(53,53),paulel(120)
c
      real*8 MMFF_ka,MMFF_ka1,MMFF_ka2
      real*8 MMFF_ka3,MMFF_ka4,MMFF_ka5
      real*8 MMFF_ka6,MMFF_ka7,MMFF_ka8
      real*8 MMFF_teta0,MMFF_teta01,MMFF_teta02
      real*8 MMFF_teta03,MMFF_teta04,MMFF_teta05
      real*8 MMFF_teta06,MMFF_teta07,MMFF_teta08
      common /merck3/ MMFF_ka(0:99,0:99,0:99),
     &                MMFF_ka1(0:99,0:99,0:99),
     &                MMFF_ka2(0:99,0:99,0:99),
     &                MMFF_ka3(0:99,0:99,0:99),
     &                MMFF_ka4(0:99,0:99,0:99),
     &                MMFF_ka5(0:99,0:99,0:99),
     &                MMFF_ka6(0:99,0:99,0:99),
     &                MMFF_ka7(0:99,0:99,0:99),
     &                MMFF_ka8(0:99,0:99,0:99),
     &                MMFF_teta0(0:99,0:99,0:99),
     &                MMFF_teta01(0:99,0:99,0:99),
     &                MMFF_teta02(0:99,0:99,0:99),
     &                MMFF_teta03(0:99,0:99,0:99),
     &                MMFF_teta04(0:99,0:99,0:99),
     &                MMFF_teta05(0:99,0:99,0:99),
     &                MMFF_teta06(0:99,0:99,0:99),
     &                MMFF_teta07(0:99,0:99,0:99),
     &                MMFF_teta08(0:99,0:99,0:99)
c
      real*8 stbn_abc,stbn_cba
      real*8 stbn_abc1,stbn_cba1
      real*8 stbn_abc2,stbn_cba2
      real*8 stbn_abc3,stbn_cba3
      real*8 stbn_abc4,stbn_cba4
      real*8 stbn_abc5,stbn_cba5
      real*8 stbn_abc6,stbn_cba6
      real*8 stbn_abc7,stbn_cba7
      real*8 stbn_abc8,stbn_cba8
      real*8 stbn_abc9,stbn_cba9
      real*8 stbn_abc10,stbn_cba10
      real*8 stbn_abc11,stbn_cba11
      real*8 defstbnd_abc,defstbnd_cba
      common /merck4/ stbn_abc(100,100,100),stbn_cba(100,100,100),
     &                stbn_abc1(100,100,100),stbn_cba1(100,100,100),
     &                stbn_abc2(100,100,100),stbn_cba2(100,100,100),
     &                stbn_abc3(100,100,100),stbn_cba3(100,100,100),
     &                stbn_abc4(100,100,100),stbn_cba4(100,100,100),
     &                stbn_abc5(100,100,100),stbn_cba5(100,100,100),
     &                stbn_abc6(100,100,100),stbn_cba6(100,100,100),
     &                stbn_abc7(100,100,100),stbn_cba7(100,100,100),
     &                stbn_abc8(100,100,100),stbn_cba8(100,100,100),
     &                stbn_abc9(100,100,100),stbn_cba9(100,100,100),
     &                stbn_abc10(100,100,100),stbn_cba10(100,100,100),
     &                stbn_abc11(100,100,100),stbn_cba11(100,100,100),
     &                defstbnd_abc(0:4,0:4,0:4),
     &                defstbnd_cba(0:4,0:4,0:4)
c
      real*8 t1_1,t2_1,t3_1
      real*8 t1_2,t2_2,t3_2
      character*16 kt_1,kt_2
      common /merck5/ t1_1(2,0:2000),t2_1(2,0:2000),t3_1(2,0:2000),
     &                t1_2(2,0:2000),t2_2(2,0:2000),t3_2(2,0:2000),
     &                kt_1(0:2000),kt_2(0:2000)
c
      real*8 G,alph,Nn
      character*1 DA
      common /merck6/ G(maxclass),alph(maxclass),Nn(maxclass),
     &                DA(maxclass)
c
      real*8 bci,bci_1
      real*8 pbci,fcadj
      common /merck7/ bci(100,100),bci_1(100,100),pbci(maxclass),
     &                fcadj(maxclass)
