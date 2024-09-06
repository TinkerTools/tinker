#
#
#  ##################################################################
#  ##                                                              ##
#  ##  compapbs.make  --  compile Tinker routines needed for APBS  ##
#  ##               (Intel Fortran for macOS Version)              ##
#  ##                                                              ##
#  ##################################################################
#
#
icc -c -O3 -no-ipo -no-prec-div -mdynamic-no-pic -w -vec-report0 pmpb.c \
 -I../apbs/include/apbs -I../apbs/include/maloc
