#
#
#  ##################################################################
#  ##                                                              ##
#  ##  compapbs.make  --  compile Tinker routines needed for APBS  ##
#  ##               (Intel Fortran for Linux Version)              ##
#  ##                                                              ##
#  ##################################################################
#
#
icc -c -O3 -axSSE3 -no-ipo -no-prec-div -w pmpb.c -I ../apbs/linux/include
