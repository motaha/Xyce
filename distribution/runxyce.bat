::  runxyce.bat
::  Sets runtime environment variables and executes Xyce 
@echo off
echo -------------------------------------------------------------------------
echo Copyright Notice
echo.
echo Copyright (c) 2002-2013, Sandia Corporation, Albuquerque, NM, USA.  
echo See the output of Xyce -license for details. 
echo -------------------------------------------------------------------------
echo.
echo.
date /t
time /t
echo.
 
"%~dp0\Xyce.exe" %*
