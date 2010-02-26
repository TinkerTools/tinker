#
#
#  ####################################################################
#  ##                                                                ##
#  ##  compgui.make  --  compile the TINKER modules needed for APBS  ##
#  ##          (Intel Fortran Compiler for Mac OSX Version)          ##
#  ##                                                                ##
#  ####################################################################
#
#
icc -c -fast -no-ipo -w -vec-report0 pmpb.c -I ../apbs/macosx/include
