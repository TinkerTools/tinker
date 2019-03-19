#
#
#  ################################################################
#  ##                                                            ##
#  ##  compgui.make  --  compile Tinker routines needed for FFE  ##
#  ##              (GNU gfortran for Linux Version)              ##
#  ##                                                            ##
#  ################################################################
#
#
#  This assumes $JAVA_HOME is set to the JAVA install directory,
#  which is /usr/lib/jvm/default-java on many Linux systems
#
#
gcc -c -O2 server.c -I $JAVA_HOME/include -I $JAVA_HOME/include/linux
