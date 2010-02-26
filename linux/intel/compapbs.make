#
#
#  ####################################################################
#  ##                                                                ##
#  ##  compgui.make  --  compile the TINKER modules needed for APBS  ##
#  ##           (Intel Fortran Compiler for Linux Version)           ##
#  ##                                                                ##
#  ####################################################################
#
#
icc -c -fast -no-ipo -w -vec-report0 pmpb.c -I ../apbs/linux/include
