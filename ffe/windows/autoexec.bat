REM
REM THIS FILE SHOULD BE APPENDED TO YOUR C:\AUTOEXEC.BAT FILE PRIOR
REM TO RUNNING FORCE FIELD EXPLORER; IF YOU DO NOT ALREADY HAVE AN
REM AUTOEXEC.BAT FILE, USE THIS ONE; PLEASE REBOOT YOUR MACHINE FOR
REM CHANGES TO AUTOEXEC.BAT TO TAKE EFFECT
REM
REM Set the three variables below to reflect local directory choices:
REM JRE is the directory containing JAVA.EXE, the Java executable
REM TINKER is the directory containing the TINKER executable files
REM FFE is the directory with FFE.JAR, Force Field Explorer's JAR file
REM
set JRE="C:\Program Files\Java\bin"
set TINKER=E:\tinker\bin
set FFEHOME=E:\tinker\jar
REM
REM You should not have to modify lines below this point
REM
set PATH=%JRE%;%tinker%;%PATH%
set CLASSPATH=%FFEHOME%\ffe.jar;%CLASSPATH%
