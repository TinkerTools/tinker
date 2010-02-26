#
#
#  #################################################################
#  ##                                                             ##
#  ##  shared.make  --  create shared object library for FFE GUI  ##
#  ##          (Intel Fortran Compiler for Linux Version)         ##
#  ##                                                             ##
#  #################################################################
#
#
#icc -c -O2 -w NativeExec.c -I /local/java/jdk1.5.0/include -I /local/java/jdk1.5.0/include/linux
#xild -shared -soname libffe.so -o libffe.so NativeExec.o
#rm NativeExec.o
#
#
#  #################################################################
#  ##                                                             ##
#  ##            (GNU gcc Compiler for Linux Version)             ##
#  ##                                                             ##
#  #################################################################
#
#
gcc -c  NativeExec.c -I /local/java/jdk1.5.0/include -I /local/java/jdk1.5.0/include/linux
gcc -shared -static NativeExec.o -o libffe.so
rm  NativeExec.o
