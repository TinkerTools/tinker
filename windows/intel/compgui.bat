@echo off
rem
rem
rem  ##############################################################
rem  ##                                                          ##
rem  ##  compgui.bat  --  compile TINKER modules needed for FFE  ##
rem  ##            (Intel Fortran for Windows Version)           ##
rem  ##                                                          ##
rem  ##############################################################
rem
rem
icl /c /O2 /w /Dwindowsintel server.c /I "C:\Program Files\Java\jdk1.6.0_22\include" /I "C:\Program Files\Java\jdk1.6.0_22\include\win32"
