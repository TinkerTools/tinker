#
#
#  ###################################################################
#  ##                                                               ##
#  ##  compgui.make  --  compile the TINKER modules needed for FFE  ##
#  ##           (Intel Fortran Compiler for Linux Version)          ##
#  ##                                                               ##
#  ###################################################################
#
#
icc -c -O3 -no-prec-div -static -w server.c -I /local/java/jdk1.5.0/include -I /local/java/jdk1.5.0/include/linux
