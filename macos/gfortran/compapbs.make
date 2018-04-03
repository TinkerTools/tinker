#
#
#  #############################################################
#  ##                                                         ##
#  ##  compapbs.make  --  compile pmpb.c using APBS headers   ##
#  ##             (GNU gfortran for MacOS Version)            ##
#  ##                                                         ##
#  #############################################################
#
#  The first include directory below (../apbs) has apbscfg.h
#  The second (/usr/local/inlcude) includes the rest following
#  completion of "make install" when buildings APBS 
#
gcc -c pmpb.c -I ../apbs-config -I /usr/local/include -I /usr/local/include/apbs

