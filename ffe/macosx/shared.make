#
#
#  #################################################################
#  ##                                                             ##
#  ##           (GNU gcc Compiler for Mac OSX Version)            ##
#  ##                                                             ##
#  #################################################################
#
#
gcc -c NativeExec.c -I /Library/Java/Home/include
gcc -dynamiclib NativeExec.o -o libffe.jnilib
rm NativeExec.o
