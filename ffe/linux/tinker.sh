#
# This is a sh script to set environment variables
# needed to use TINKER with the Force Field Explorer GUI
#
# The values of the JRE, TINKER and FFEHOME directories
# below must be modified to reflect your local choices:
# JRE is the directory containing the Java executable
# TINKER is the directory containing TINKER executables
# FFEHOME is the directory with ffe.jar and libffe.so
#

#
# Set location of JRE, TINKER and Force Field Explorer
#
JRE="/local/java/jdk1.5.0/jre"
TINKER="/user/ponder/tinker"
FFEHOME="/user/ponder/tinker/jar"
export JRE
export TINKER
export FFEHOME

#
# Append Java and TINKER binaries to default PATH
#
if [ -z "${PATH}" ]
then
   PATH="$TINKER/bin:$JRE/bin"
else
   PATH="$TINKER/bin:$JRE/bin:$PATH"
fi
export PATH

#
# Append Force Field Explorer JAR file to CLASSPATH
#
if [ -z "${CLASSPATH}" ]
then
   CLASSPATH="$FFEHOME/ffe.jar"
else
   CLASSPATH="$FFEHOME/ffe.jar:$CLASSPATH"
fi
export CLASSPATH

#
# Append JRE Shared Object Libraries to Library Path
#
if [ -z "${LD_LIBRARY_PATH}" ]
then
   LD_LIBRARY_PATH="$JRE/lib/i386/client:$JRE/lib/i386:$TINKER/bin:$FFEHOME"
else
   LD_LIBRARY_PATH="$JRE/lib/i386/client:$JRE/lib/i386:$TINKER/bin:${FFEHOME}:$LD_LIBRARY_PATH"
fi
export LD_LIBRARY_PATH

#
# Create an alias to execute Force Field Explorer
#
alias ffe="java -server -mx256M -Dtinker.dir=$TINKER -Djava.library.path=$FFEHOME -jar $FFEHOME/ffe.jar"
