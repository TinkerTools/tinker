@echo off
set FFEHOME=E:\tinker\jar
java -server -mx256M -Djava.library.path=%FFEHOME% -jar %FFEHOME%\ffe.jar
