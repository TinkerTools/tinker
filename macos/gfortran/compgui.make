#
#
#  ################################################################
#  ##                                                            ##
#  ##  compgui.make  --  compile Tinker routines needed for FFE  ##
#  ##              (GNU gfortran for macOS Version)              ##
#  ##                                                            ##
#  ################################################################
#
#
#  This assumes $JAVA_HOME is set to the JAVA install directory,
#  found under /Library/Java/JavaVirtualMachines on macOS systems
#
#
gcc -c -O2 server.c -I $JAVA_HOME/include -I $JAVA_HOME/include/linux
