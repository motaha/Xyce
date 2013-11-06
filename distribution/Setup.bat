::  Setup.bat
:: 

::  hide batch file commands
@echo off

::  clear screen
cls


::  check for the default cygwin setup 
set CYGDEF=%SYSTEMDRIVE%\cygwin


::  check for cygwin dll loaded in memory
set running=XYCLEAN

:: DEBUG:  2010.12.06 hacked to enable regression runs on Win7
:: DEBUG:  DO NOT DISTRIBUTE!!!!!!!!! Use "make package" to create native installer
:: for /f "tokens=1" %%y in ('tasklist /M "cygwin*" /NH') do set running=%%y

::  clear screen
cls


::  ignore INFO message when no cygwin dll is loaded
if %running% == INFO: set running=XYCLEAN

:: DEBUG:  2010.12.06 hacked to enable regression runs on Win7
:: DEBUG:  DO NOT DISTRIBUTE!!!!!!!!! Use "make package" to create native installer
set running=DEBUG


if %running% == XYCLEAN ( 

  ::  change to installer directory and execute script with local cygwin
  cd %~dp0
  set CYGLOCAL=1
  ".\files\cygbin\bash.exe" -e ".\files\install4win.sh" %* 2>NUL 

) else if exist %CYGDEF% (
  ::  change to installer directory and execute script with system cygwin
  cd %~dp0
  set CYGLOCAL=0
  "%CYGDEF%\bin\bash.exe" -e ".\files\install4win.sh" %* 2>NUL 

) else (

  :: alert the user to the problem
  echo Cygwin is not installed in the default location:  %CYGDEF%
  echo.
  echo and the following running Cygwin applications conflict with Setup: 
  echo. 
  for /f "tokens=1" %%y in ('tasklist /M "cygwin*" /NH') do echo 	%%y
  echo. 
  echo Please close the programs and try again. 
  echo.

  :: set errorlevel
  echo a | find "b" >NUL
)



::  handle errors
if errorlevel 1 (
  echo.
  echo Setup did not complete successfully.  
  echo.
  echo Consult the Xyce User's Guide for help with installation. 
  echo.
  echo.
) else (
  echo.
  echo Xyce installation completed successfully.
  echo.
)


::  press key or close cmd line window
pause

