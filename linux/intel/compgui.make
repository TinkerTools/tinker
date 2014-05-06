#
#
#  ################################################################
#  ##                                                            ##
#  ##  compgui.make  --  compile TINKER routines needed for FFE  ##
#  ##              (Intel Fortran for Linux Version)             ##
#  ##                                                            ##
#  ################################################################
#
#
icc -c -O3 -no-prec-div -static -w server.c -I /local/java/jdk1.5.0/include -I /local/java/jdk1.5.0/include/linux
