#
# This is a csh/tcsh script to set environment variables
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
setenv JRE "/local/java/jdk1.5.0/jre"
setenv TINKER "/user/ponder/tinker"
setenv FFEHOME "/user/ponder/tinker/jar"
setenv LD_ASSUME_KERNEL 2.2.5

#
# Append Java and TINKER binaries to default PATH
#
if !($?PATH) then
   setenv PATH "$TINKER/bin:$JRE/bin"
else
   setenv PATH "$TINKER/bin:$JRE/bin:$PATH"
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
# Append JRE Shared Object Libraries to Library Path
#
if !($?LD_LIBRARY_PATH) then
   setenv LD_LIBRARY_PATH "$JRE/lib/i386/client:$JRE/lib/i386:$TINKER/bin:$FFEHOME"
else
   setenv LD_LIBRARY_PATH "$JRE/lib/i386/client:$JRE/lib/i386:$TINKER/bin:${FFEHOME}:$LD_LIBRARY_PATH"
endif

#
# Create an alias to execute Force Field Explorer
#
alias ffe "java -server -mx256M -Dtinker.dir=$TINKER -Djava.library.path=$FFEHOME -jar $FFEHOME/ffe.jar"
