@echo off
rem
rem
rem  ###########################################################
rem  ##                                                       ##
rem  ##  dll.bat  --  build ffe.dll for Force Field Explorer  ##
rem  ##     (Intel Fortran Compiler for Windows Version)      ##
rem  ##                                                       ##
rem  ###########################################################
rem
rem
icl /LD /w NativeExec.c /o ffe.dll /I "c:\Program Files\Java\jdk1.5.0_01\include" /I "c:\Program Files\Java\jdk1.5.0_01\include\win32"
del NativeExec.obj ffe.lib ffe.exp
rem
rem
rem  ###########################################################
rem  ##                                                       ##
rem  ##        (GNU gcc Compiler for Windows Version)         ##
rem  ##                                                       ##
rem  ###########################################################
rem
rem
rem gcc -c NativeExec.c -o NativeExec.o -I "c:\Program Files\Java\jdk1.5.0_01\include" -I "c:\Program Files\Java\jdk1.5.0_01\include\win32"
rem dlltool NativeExec.o -A --export-all-symbols --output-def ffe.def
rem dllwrap NativeExec.o -static --def ffe.def -o ffe.dll
