#!/bin/sh

#   install.sh
#   Extracts tarball into user selected directory and configures shell script for
#   running Xyce from command line.  


echo "-----------------------------------------------------------------------------"
echo " Copyright Notice"
echo 
echo " Copyright (c) 2002-2013, Sandia Corporation, Albuquerque, NM, USA.  "
echo " See the output of Xyce -license for details. "
echo "-----------------------------------------------------------------------------"
echo
echo
echo "Preparing to install Xyce..."


#   setup path to needed progs
PATH="$PATH:/bin:/usr/bin:/usr/local/bin:/sbin:/usr/sbin"


# set envvars 
VER=5.3
#ORIG_PWD=`pwd`

ORIG_PWD=`dirname $0`
echo "ORIG_PWD is $ORIG_PWD"
if [ "x$ORIG_PWD" = "x." ]
then
   ORIG_PWD=`pwd`
fi

# trap errors and ctrl-c 
alert_user ()
{
  echo
  echo "Setup interrupted..."
  $CNOK && echo "Installation aborted:  alert_user() called" >> "$CONFIG_NOTES"
  exit 11;
}
trap "alert_user" 1 2 3 4 6 8 10 12 13 15


#   make sure we have the install tarball
if [ ! -f "$ORIG_PWD/xtl.tar.gz" ]; then
    echo "ERROR:  There are missing install files.  Please unpack a"
    echo "ERROR:  valid install tarball and run install.sh again."
    echo
    echo "ERROR:  You must run this install script from the"
    echo "ERROR:  directory where it is unpacked."
    echo
    echo "ERROR:  Installation aborted."
    exit 1;
fi


#   record configuration options
CONFIG_NOTES="$ORIG_PWD/install.log.$$"
echo "Xyce Installation Notes" > "$CONFIG_NOTES"


#   what plat is this?  (in case we have a mega install bundle :)
OS_NAME=`uname -s 2> /dev/null | tr "[:upper:]" "[:lower:]" 2> /dev/null`
echo "platform:  $OS_NAME" >> "$CONFIG_NOTES"


# set default prefix
INSTALLBASE="$HOME/Xyce" 
DONOTPROMPT=0

if [ $# = 0 ]
then
  echo "No command line args given"
else
  INSTALLBASE=$1
  DONOTPROMPT=1
fi
    
#   determine where to install Xyce    
INVALID="INVALID"
while [ "$INVALID" = "INVALID" ]
do
    if [ $DONOTPROMPT = 0 ]
    then
      # a ? means the last entry was bogus
      echo "Where should Xyce be installed?  [ $INSTALLBASE ]  "
      read ANS  
    
      #   user typed something so expand and check for validity
      if [ "x$ANS" != "x" ]; then
          INSTALLBASE=`eval echo $ANS`
      fi
    fi
   
    #   throw out bogus entries
    WALKER="DONE"
    case "$INSTALLBASE" in
        *\?* | *\** | *\$* )
            echo "WARNING:  $INSTALLBASE is not a valid path name."
	    echo
        ;;
        
	*\ * )
	    echo "WARNING:  Path must not contain spaces."
	;;

        * )
            WALKER="$INSTALLBASE"
        ;;
    esac
  
    #   make sure path is writeable

    while [ "$WALKER" != "DONE" ]
    do
        if [ -d "$WALKER" ]; then
            if [ ! -w "$WALKER" ]; then
                #   cannot write to this path so alert user and try again
                echo "WARNING:  $INSTALLBASE is not writeable."
                echo "WARNING:  Please select a new location."
                INSTALLBASE="?"
            else
                #   path is writeable and we can proceed to copy files
                INVALID="ok"
            fi
            #   reached root (/) so stop climbing
            WALKER="DONE" 
        else
            #   did not find an existing file so step up the hierarchy
            #   removing the trailing dirname and preceding slash
            WALKER=`echo $WALKER | sed 's/[/]$//' | sed 's/[^/]*$//'`
            
            #   found null walker so create the path in (relative to) pwd 
            if [ "x$WALKER" = "x" ]; then
                INVALID="ok"
                WALKER="DONE"
            fi
        fi
    done
    
    #   dbl check for accidentally writing over another installation
    #   skip this test if the path is not writeable
    if [ "$INVALID" = "ok" -a -d "$INSTALLBASE" ]; then
        while [ "$ANS" != "y" -a "$ANS" != "Y" -a "$ANS" != "n" -a "$ANS" != "N" ]
        do
            if [ $DONOTPROMPT = 0 ]
            then
              echo "WARNING:  $INSTALLBASE already exists.  Some files may be overwritten."
              echo "WARNING:  Do you wish to continue?  [y/n]"
              read ANS
            else
              ANS="Y"
            fi
        done
        if [ "$ANS" = "N" -o "$ANS" = "n" ]; then
            INVALID="INVALID"
            INSTALLBASE="?"
        fi
    fi 
    if [ $INVALID = "INVALID" -a $DONOTPROMPT = 1 ]
    then
       exit 1
    fi  
done

echo "Xyce will be installed in $INSTALLBASE"

echo "Decompressing archive and copying files..."
mkdir -p "$INSTALLBASE"
if [ "$?" -ne "0" ]; then
    echo "ERROR:  There was a problem creating $INSTALLBASE.  Exiting." 
    exit 1
fi

cd "$INSTALLBASE"

INSTALLBASE=`pwd`

cp -f "$ORIG_PWD/xtl.tar.gz" "$INSTALLBASE"

if [ "$?" -ne "0" ]; then
    echo "ERROR:  There was a problem copying files to the installation directory.  Exiting." 
    exit 1
fi

gzip -d "$INSTALLBASE/xtl.tar.gz"
if [ "$?" -ne "0" ]; then 
    echo "ERROR:  There was a problem decompressing the archive.  Exiting."
    exit 1
fi

tar xf "$INSTALLBASE/xtl.tar" > /dev/null 2>&1
if [ "$?" -ne "0" ]; then
    echo "ERROR:  There was a problem extracting files from the archive.  Exiting." 
    exit 1
fi


#   clean up (note that if user tries to install in $ORIG_PWD the tarball
#   gets wiped.  user can recover by unpacking installer again
rm -f "$INSTALLBASE/xtl.tar"

echo "installing Xyce here:   $INSTALLBASE" >> "$CONFIG_NOTES"


#   create the executable script (runxyce for serial and xmpirun for parallel)
echo "Configuring Xyce for this computer..."

case $0 in
    *MPI* | *mpi* ) 
        XYCERUN="xmpirun" 
    ;;
    
    * )             
        XYCERUN="runxyce" 
    ;;
esac


XYCESTART="$INSTALLBASE/bin/$XYCERUN"


echo "#!/bin/sh" > "$XYCESTART"
echo " " >> "$XYCESTART"

echo "echo \"-----------------------------------------------------------------------------\"" >> "$XYCESTART"
echo "echo \" Copyright Notice\"" >> "$XYCESTART"
echo "echo" >> "$XYCESTART"
echo "echo \" Copyright (c) 2002-2013, Sandia Corporation, Albuquerque, NM, USA.  \"" >> "$XYCESTART"
echo "echo \" See the output of Xyce -license for details. \"" >> "$XYCESTART"
echo "echo \"-----------------------------------------------------------------------------\"" >> "$XYCESTART"
echo "echo" >> "$XYCESTART"
echo "echo" >> "$XYCESTART"
echo "date" >> "$XYCESTART"
echo "echo" >> "$XYCESTART"
echo " " >> "$XYCESTART"
echo "LOCDIR=\"\$( cd \"\$( dirname \"\$0\" )\"/.. && pwd )\"" >> $XYCESTART


# set lib env vars
case $0 in 

  *Linux*|*RHEL6* )
    echo "export LD_LIBRARY_PATH=\"\$LOCDIR/lib:\$LD_LIBRARY_PATH\"" >> "$XYCESTART"
  ;;

  *OSX* )
    echo "export DYLD_LIBRARY_PATH=\"\$LOCDIR/lib:\$DYLD_LIBRARY_PATH\"" >> "$XYCESTART"
  ;;

  * )
    # nothing to do
  ;;

esac


#   wrap Xyce executable in proper mpirun call if necessary
#   FIXME:  need to detect the current mpi env and configure
if [ "$XYCERUN" = "xmpirun" ]; then
    #   mpirun script depends on platform and mpi implementation
    case $0 in
	*TLCC* | *NWCC* | *Tbird*  )
            echo "# Use system modules to set environment and path.  See module help pages." >> "$XYCESTART"
            echo "# Note that executing this script may change the current environment settings." >> "$XYCESTART"
            ## module loads do not return exit codes so alert users to process
            echo "echo \"Load module for Xyce $VER using:  \"" >> "$XYCESTART"

            ## skipping version of modules due to limited access to module setup
            echo "echo \"               module load xyce\"" >> "$XYCESTART"

            ## mpirun / mpiexec are symlinks to orterun in Open MPI >= 1.2.5
	    MPISCRIPT="mpirun"
        ;;
        
        *OPENMPI* )
            # Using default prefix for openmpi (/usr/local)
            echo "# Uncomment and set PATH to the openmpi binaries here:" >> "$XYCESTART"
	    echo "#  PATH=\"/opt/openmpi/bin:\$PATH\"" >> "$XYCESTART"
	    MPISCRIPT="mpirun"
        ;;

        *MPICH* )
            # Using default prefix for mpich/mpich2 (/usr/local/mpich)
            echo "# Uncomment and set PATH to the mpich binaries here:" >> "$XYCESTART"
	    echo "#  PATH=\"/usr/local/mpich/bin:\$PATH\"" >> "$XYCESTART"
	    MPISCRIPT="mpirun"
        ;;

        * )
            echo "# Set MPISCRIPT to the full path to mpirun here:" >> "$XYCESTART"
	    MPISCRIPT="mpirun"
        ;;
    esac

    # write the var to file
    echo "MPISCRIPT=\"$MPISCRIPT\"" >> "$XYCESTART"
    echo " " >> "$XYCESTART"
    echo " " >> "$XYCESTART"
    echo " " >> "$XYCESTART"
    echo " " >> "$XYCESTART"
    

    # at install time, alert user if mpi not found
    [ ! -x $MPISCRIPT ] && echo "The MPI environment could not be determined.  This may require additional configuration prior to running Xyce with xmpirun." && echo "$MPISCRIPT not found/executable" >> $CONFIG_NOTES 


    # at run time, check for valid mpirun/mpiexec/etc
    echo "which \$MPISCRIPT >/dev/null 2>&1" >> $XYCESTART
    echo "[ \$? -ne 0 ] && echo \"Xyce failed to run because the MPI environment could not be determined.  Please edit xmpirun, or update your path to include mpirun/mpiexec/etc.  Consult the Xyce Reference Guide for assistance.\" && exit 1" >> "$XYCESTART"

    
    # at run time,  preprocess $MPISCRIPT and Xyce options                          
    echo "until [ \$# -eq 0 ]" >> "$XYCESTART"
    echo "do" >> "$XYCESTART"
    echo "    case \$1 in" >> "$XYCESTART"

    echo "        -np )" >> "$XYCESTART"
    echo "            # push the -np option to the mpi run script" >> "$XYCESTART" 
    echo "            MPIRUNARGS=\"\$MPIRUNARGS \$1 \$2\"" >> "$XYCESTART"    
    echo "            shift" >> "$XYCESTART"
    echo "            shift" >> "$XYCESTART"
    echo "        ;;" >> "$XYCESTART"

    echo "        -v )" >> "$XYCESTART"
    echo "            # short circuit to keep -v from leaking to mpirun" >> "$XYCESTART" 
    echo "            \$MPISCRIPT -np 1 \"\$LOCDIR/bin/Xyce\" -v" >> "$XYCESTART"
    echo "            exit 0" >> "$XYCESTART"
    echo "        ;;" >> "$XYCESTART"

    echo "        -h )" >> "$XYCESTART"
    echo "            echo \"Please consult the $MPISCRIPT documentation for more options.\"" >> "$XYCESTART"
    echo "            echo" >> "$XYCESTART"
    echo "            echo" >> "$XYCESTART"
    echo "            \$MPISCRIPT -np 1 \"\$LOCDIR/bin/Xyce\" -h" >> "$XYCESTART"
    echo "            echo" >> "$XYCESTART"  
    echo "            echo" >> "$XYCESTART"  
    echo "            echo \"Usage:  xmpirun [$MPISCRIPT options] -np <processors> [Xyce options] <netlist filename>\"" >> "$XYCESTART"   
    echo "            exit 0" >> "$XYCESTART"
    echo "        ;;" >> "$XYCESTART"

    echo "        * )" >> "$XYCESTART"
    echo "            XYCEARGS=\"\$XYCEARGS \$1\"" >> "$XYCESTART"
    echo "            shift" >> "$XYCESTART"    
    echo "        ;;" >> "$XYCESTART"

    echo "    esac" >> "$XYCESTART"

    echo "done" >> "$XYCESTART"
    echo "exec \$MPISCRIPT \$MPIRUNARGS \"\$LOCDIR/bin/Xyce\" \$XYCEARGS" >> "$XYCESTART"
else
    echo "exec \"\$LOCDIR/bin/Xyce\" \$*" >> "$XYCESTART"    
fi


chmod 755 "$XYCESTART"
echo "created executable script:  $XYCESTART" >>  "$CONFIG_NOTES"


#   create symlink to Xyce launher
XYCESYM=$HOME/$XYCERUN

#Comment this out: 
#   preserve old executable in case user had some real program there
#if [ -f $XYCESYM ]; then
#    echo "copied existing $XYCESYM to $XYCESYM.save.$$" >> "$CONFIG_NOTES"
#    mv $XYCESYM $XYCESYM.save.$$
#fi
# and do this instead.  These save files are starting to fill up home 
# directories!  This script is deprecated anyway and now only used in testing,
# so who cares about saving the junk.
rm -f $XYCESYM

ln -sf "$XYCESTART" $XYCESYM
echo "created symlink here:  $XYCESYM" >> "$CONFIG_NOTES"


cd "$ORIG_PWD"
echo "$CONFIG_NOTES contains configuration notes."
echo
echo "Xyce installation completed successfully."
echo 

exit 0
