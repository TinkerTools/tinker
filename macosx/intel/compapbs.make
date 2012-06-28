#
#
#  #####################################################################
#  ##                                                                 ##
#  ##  compapbs.make  --  compile the TINKER modules needed for APBS  ##
#  ##           (Intel Fortran Compiler for Mac OSX Version)          ##
#  ##                                                                 ##
#  #####################################################################
#
#
icc -c -O3 -no-ipo -no-prec-div -mdynamic-no-pic -w -vec-report0 pmpb.c -I ../apbs/macosx/include
