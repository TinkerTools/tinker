#
# This is a sh script to set environment variables
# needed to use TINKER with the Force Field Explorer GUI
#
# The values of the TINKER and FFEHOME directories
# below must be modified to reflect your local choices:
# TINKER is the directory containing TINKER executables
# FFEHOME is the directory with ffe.jar and libffe.so
#

#
# Set location of JRE, TINKER and Force Field Explorer
#
TINKER="/user/ponder/tinker"
FFEHOME="/user/ponder/tinker/jar"
export TINKER
export FFEHOME

#
# Append the TINKER binaries to default PATH
#
if [ -z "${PATH}" ]
then
   PATH="$TINKER/bin"
else
   PATH="$TINKER/bin:$PATH"
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
# Create an alias to execute Force Field Explorer
#
alias ffe="java -server -mx256M -Dtinker.dir=$TINKER -Djava.library.path=$FFEHOME -jar $FFEHOME/ffe.jar"
