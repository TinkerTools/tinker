#
#
#  #####################################################################
#  ##                                                                 ##
#  ##  compapbs.make  --  compile the TINKER modules needed for APBS  ##
#  ##            (Intel Fortran Compiler for Linux Version)           ##
#  ##                                                                 ##
#  #####################################################################
#
#
icc -c -O3 -axSSE3 -no-ipo -no-prec-div -w pmpb.c -I ../apbs/linux/include
