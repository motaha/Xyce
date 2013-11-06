#!/bin/bash

# install4win.sh
# Cygwin derived shell script for installing and configuring Xyce


# look for cygwin environment settings using shell builtins
case $CYGLOCAL in
  1 )  
    CYGWINDLLVER=1.5.12-1
    PATH="`pwd`/files/cygbin"
    ;;
  0 )
    CYGWINDLLVER="SYSTEM"
    export PATH=/bin:/usr/bin:/usr/local/bin:/sbin
    ;;
  * )
    # error exit now to prevent dll conflict
     exit 1
    ;;
esac


# set envvars 
VER=5.3
ORIG_PWD=`pwd`  
PLAT=`uname -a` 
CONFIG_NOTES="$ORIG_PWD/configuration notes $$.txt"


# trap errors and ctrl-c 
alert_user () 
{
  echo
  echo "Setup interrupted..."
  $CNOK && echo "Installation aborted:  alert_user() called" >> "$CONFIG_NOTES"
  exit 11;
}
trap "alert_user" 1 2 3 4 6 8 10 12 13 15



# clean up
echo
echo



# test for writing configuration notes to top level dir 
if [ -w "$ORIG_PWD" ]; then
  CNOK="true"
else
  echo
  echo "$ORIG_PWD is not writeable.  Installation notes will not be saved."
  CNOK="false"
fi



# open installer log file
$CNOK && echo "Xyce Configuration Notes" > "$CONFIG_NOTES"
$CNOK && echo "platform:  $PLAT" >> "$CONFIG_NOTES"



# display copyright message
echo
echo "-------------------------------------------------------------------------"
echo "Copyright Notice"
echo 
echo "Copyright (c) 2002-2013, Sandia Corporation, Albuquerque, NM, USA. "
echo " See the output of Xyce -license for details. "
echo "-------------------------------------------------------------------------"



# begin installation
echo
echo "Preparing to install Xyce..."



# find installation files 
if [ ! -r "$ORIG_PWD/files/xtl.tar.gz" ]; then
  $CNOK && echo "Installation aborted:  xtl not readable " >> "$CONFIG_NOTES"
  echo 
  echo "Setup cannot continue.  There are missing or unreadable files."
  exit 1;
fi



# check for administrator read/write privileges
if [ "x$SYSTEMDRIVE" != "x" -a -w "$SYSTEMROOT" ]; then
  INSTALLBASE=`cygpath -ua "$SYSTEMDRIVE/Xyce-$VER"`
else
  if [ -d "$USERPROFILE" -a -w "$USERPROFILE" ]; then
    INSTALLBASE=`cygpath -ua "$USERPROFILE\Xyce-$VER"`
  else
    INSTALLBASE=""; 
  fi
fi
DONOTPROMPT=0

# allow command line argument to override default path and disable prompting
if [ $# != 0 ]
then
  INSTALLBASE=$1
  DONOTPROMPT=1
fi

# determine where to install Xyce    
if [ "x$INSTALLBASE" != "x" ]; then
  SHOWDEFAULT="true"
fi



# get user input
INVALID="INVALID"



# display path in Windows format to the user
if [ "x$INSTALLBASE" != "x" ]; then
  WINSTALLBASE=`cygpath -wla "$INSTALLBASE"`
fi



# request user input
while [ "$INVALID" = "INVALID" ]
do

  if [ $DONOTPROMPT = 0 ]
  then   
    # clean up
    echo

    # show input prompt
    if [ "x$SHOWDEFAULT" = "xtrue" ]; then
      SHOWDEFAULT="false"
      echo "Please type the name of the folder where Xyce $VER should be installed, or press ENTER to select the default location."
      echo
      echo -n "[ $WINSTALLBASE ]  " 
    else
      echo "Please type the name of the folder where Xyce $VER should be installed."
      echo
      echo -n "?  "
    fi
  
    # get user input
    read -r "ANS"  

    # ignore blank entry
    if [ "x$ANS" != "x" ]; then
      INSTALLBASE="$ANS"
  
      # convert paths internally
      INSTALLBASE=`cygpath -ua "$INSTALLBASE"`
      WINSTALLBASE=`cygpath -wla "$INSTALLBASE"`
    fi
  fi

  # throw out bogus entries 
  case "$INSTALLBASE" in
    *[\:\*\?\"\<\>]* )
      echo
      WINSTALLBASE=`cygpath -wla "$INSTALLBASE"`
      echo "$WINSTALLBASE is not a valid folder name."
      INSTALLBASE=""
    ;;
  esac

  # check for legal path
  if [ "x$INSTALLBASE" != "x" ]; then
    INVALID="ok"
  fi

  # try to copy files
  if [ "$INVALID" = "ok" ]; then

    # make sure we aren't accidentally writing over another installation
    if [ -d "$INSTALLBASE" ]; then
      if [ $DONOTPROMPT = 0 ]
      then 
        echo
        echo "The $WINSTALLBASE folder already exists.  Some files may be overwritten if you proceed."  
        echo
        echo -n "Do you wish to continue? [y/N]  "

        # drop all answers except positive confirmation
        read "ANS"
      else
         # assume if the user gave a command line arg, that he meant us to
         # do this
         ANS="Y"
      fi
      case "$ANS" in
      [yY] | [yY][eE] | [yY][eE][sS] )
      ;;
      * ) 
        INSTALLBASE=""
        INVALID="INVALID"
      ;;
      esac
      # make sure existing partial path is writeable
    else
      PATHCHECK="$INSTALLBASE"
      
      while [ "x$PATHCHECK" != "x." ]
      do
        # strip off trailing directory name
        PATHCHECK=`dirname "$PATHCHECK"`

        if [ -d "$PATHCHECK" ]; then     
          # alert user if remaining path is not writeable
          if [ ! -w "$PATHCHECK" ]; then
            INSTALLBASE=""
            INVALID="INVALID"
            echo
            echo "The folder $WINSTALLBASE could not be created.  You may not have write privileges, or the location may contain Shortcuts."
          fi          

          # path found so exit loop and continue installing
          PATHCHECK="."
        fi
      done
    fi 

    # create directory
    if [ "$INVALID" = "ok" ]; then
      mkdir -p "$INSTALLBASE"
      
      # final check against existing top level dir w/o r/w perms
      if [ ! -w "$INSTALLBASE" ]; then
        INSTALLBASE=""
        INVALID="true"
        echo
        echo "Unable to save files in $WINSTALLBASE.  Access is denied.  Please select an alternate location."
      fi
   fi

  fi

  # else request location from the user again

done



echo
echo
echo "Xyce $VER will be installed in $WINSTALLBASE"
echo "     Decompressing archive and copying files..."


# avoid long filename bug in the runxyce script
WINSTALLBASE=`cygpath -wa "$INSTALLBASE"`

$CNOK && echo "installation directory:   $WINSTALLBASE" >> "$CONFIG_NOTES"

cd "$INSTALLBASE"

cp -f "$ORIG_PWD/files/xtl.tar.gz" "$INSTALLBASE"
if [ $? -ne 0 ]; then
  $CNOK && echo "problem writing to installbase" >> $CONFIG_NOTES
  exit 1
fi

gzip -d "$INSTALLBASE/xtl.tar.gz"
if [ $? -ne 0 ]; then 
  $CNOK && echo "problem decompressing the archive" >> $CONFIG_NOTES
  exit 1
fi

tar xf "$INSTALLBASE/xtl.tar" 
if [ $? -ne 0 ]; then
  $CNOK && "problem extracting files from the archive" >> $CONFIG_NOTES 
  exit 1
fi

echo "     Removing temporary files..." 
rm -f "$INSTALLBASE/xtl.tar"



# create the runxyce.bat batch file
echo "     Configuring Xyce for this computer..."
XYCESTART="$INSTALLBASE/bin/runxyce.bat"



echo "::  runxyce.bat" > "$XYCESTART"
echo "::  Sets runtime environment variables and executes Xyce " >> "$XYCESTART"
echo "@echo off" >> "$XYCESTART"
echo "cls" >> "$XYCESTART"
echo "echo -------------------------------------------------------------------------" >> "$XYCESTART"
echo "echo Copyright Notice" >> "$XYCESTART"
echo "echo." >> "$XYCESTART"
echo "echo Copyright (c) 2002-2013, Sandia Corporation, Albuquerque, NM, USA.  " >> "$XYCESTART"
echo "echo See the output of Xyce -license for details." >> "$XYCESTART"
echo "echo -------------------------------------------------------------------------" >> "$XYCESTART"
echo "echo." >> "$XYCESTART"
echo "echo." >> "$XYCESTART"
echo "date /t" >> "$XYCESTART"
echo "time /t" >> "$XYCESTART"
echo "echo." >> "$XYCESTART"
## this line properly handles quoted input, but masks Xyce's error handling
## echo "IF [%1]==[] goto error" >> "$XYCESTART"
echo " " >> "$XYCESTART"
## echo "\"$WINSTALLBASE\\bin\\Xyce\" %*" >> "$XYCESTART"
echo "\"%~dp0\\Xyce.exe\" %*" >> "$XYCESTART"
echo "goto done" >> "$XYCESTART"
echo ":ERROR" >> "$XYCESTART"
## echo "\"$WINSTALLBASE\\bin\\Xyce\" -h" >> "$XYCESTART"
echo "\"%~dp0\\Xyce.exe\" -h" >> "$XYCESTART"
echo "echo." >> "$XYCESTART"
echo "echo runxyce.bat should be executed from the command line. (cmd.exe)" >> "$XYCESTART"
echo "echo." >> "$XYCESTART"
echo "echo Consult the Xyce $VER User's Guide for assistance." >> "$XYCESTART"
echo "pause" >> "$XYCESTART"
echo ":DONE" >> "$XYCESTART"
echo "exit %ERRORLEVEL%" >> "$XYCESTART"




$CNOK && echo "created batch file (cygpath):  $XYCESTART" >>  "$CONFIG_NOTES"


echo "     Setting permissions..."
chmod +x "$XYCESTART"



# echo "     Converting file format..."
# unix2dos -D "$XYCESTART"



echo "     Confirming setup..." `"$INSTALLBASE/bin/Xyce" -v` 
echo
echo



exit 0
