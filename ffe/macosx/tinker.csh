#
# ######################################################
# ##                                                  ##
# ##  TINKER/FFE on MacOSX/csh Initialization Script  ##
# ##                                                  ##
# ######################################################
#
# This is a csh/tcsh script to set environment variables
# needed to use TINKER with the Force Field Explorer GUI
#
# The values of the TINKER and FFEHOME directories
# below must be modified to reflect your local choices:
# TINKER is the directory containing TINKER executables
# FFEHOME is the directory with ffe.jar and libffe.jnilib
#

#
# Set location of JRE, TINKER and Force Field Explorer
#
setenv TINKER "/Users/ponder/tinker"
setenv FFEHOME "/Users/ponder/tinker/jar"

#
# Append the TINKER binaries to default PATH
#
if !($?PATH) then
   setenv PATH "$TINKER/bin"
else
   setenv PATH "$TINKER/bin:$PATH"
endif

#
# Append Force Field Explorer JAR file to CLASSPATH
#
if !($?CLASSPATH) then
   setenv CLASSPATH "$FFEHOME/ffe.jar"
else
   setenv CLASSPATH "$FFEHOME/ffe.jar:$CLASSPATH"
endif

#
# Create an alias to execute Force Field Explorer
#
alias ffe "java -server -mx256M -Dtinker.dir=$TINKER -Djava.library.path=$FFEHOME -jar $FFEHOME/ffe.jar"
