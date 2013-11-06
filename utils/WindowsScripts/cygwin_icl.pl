#!/usr/bin/perl
#
# usage:
#
# cygwin_icl <compiler> <compiler args>
#
# This script was originally developed to get Intel compilers to configure
# Trilinos.  The major problem addressed was the format required by Intel 
# compilers to name the object file, if the name was different than the source
# file.  This perl script parses the command line and converts only the 
# necessary arguments, passing the rest through.  In fact, the compiler is
# an argument to this script, so it can be used elsewhere if needed.
#
# Author : Heidi Thornquist
# Date : April 7, 2003
# Revised : March 31, 2005
# Note:  This script was updated to be compatable with Intel v8 compilers.
#
# Author: Rich Schiek
# Revised : Feb. 21, 2006 
# If an absolute path is passed in for an include directory as in
# -I/home the compiler can't find this as /home doesn't exist
# outside of cygwin.  Simple fix is the change -I/home to
# -Ic:\\cygwin\\home.  This assumes cygwin is installed on c:
#
# Author: Rich Schiek
# Revised: Feb 16, 2007
# If the name of the source file comes after a linker option
# on the command line, then the compiler (ifort) can't find
# the source.  So, we'll move the source to the first arg
# after the compiler.

$compiler = $ARGV[0] . " ";
$output = "";
$i = 1;
$numinputs = @ARGV;
while ($i < $numinputs) 
{
  # First things first, fix the -o problem!
  if (index($ARGV[$i],"-o")==0 && index($ARGV[$i+1],"exe")<0) 
  {
    $output = $output . "-Fo" . $ARGV[$i+1] . " ";
    $i++;
  } 
  elsif (index($ARGV[$i],"\\")>=0) 
  {
    # Whoops, we have a little problem with windows paths :)
    $wherebeg=0; $where=0;  
    $newfilename = "";
    while ($where >=0) 
    {
      $where = index($ARGV[$i],"\\", $wherebeg );
      if ($where < 0) 
      {
        $part = substr($ARGV[$i], $wherebeg, 100);
        $newfilename = $newfilename . $part;
      } 
      else 
      {
        $part = substr($ARGV[$i], $wherebeg, $where-$wherebeg);
        $newfilename = $newfilename . $part . "\\\\";
      }
      $wherebeg = $where + 1;
    }
    $output = $output . $newfilename . " ";
  } 
  elsif (index($ARGV[$i],"-l")==0) 
  {
    # Fix the -l problem
    $newlib = "";
    $arg_length = length($ARGV[$i]);
    if ($arg_length > 2) 
    {
      # Check rest of argument, library is tacked on (ex. -lteuchos)
      $templib = substr($ARGV[$i], 2, 100);
      if (index($templib,".lib")>=0) 
      {
        $newlib = $templib;
      } 
      else 
      {
        $newlib = "lib" . $templib . ".lib";  # lib<library_name>.lib
      }
    } 
    else 
    {
      $newlib = $ARGV[$i+1];
      $i++; 
    }
    $output = $output . " " . $newlib . " ";
  } 
  elsif (index($ARGV[$i],"-L")==0) 
  {
    # Fix the -L problem
    $newlibpath = "";
    $pathpart = substr($ARGV[$i], 2, 100);
    # Check to see if this is a cygwin path (ex. /cygdrive/c/Trilinos/...)
    $wherebeg=index($ARGV[$i],"/cygdrive");
    if ($wherebeg >= 0) 
    {
      # Grab the name of the disk, which is expected to be the next directory
      # after cygdrive (ex. /cygdrive/c/Trilinos/..., "c" is the disk name)
      $part = "";
      $where=index($ARGV[$i],"/",$wherebeg+10);
      $newlibpath=$newlibpath . substr($ARGV[$i], $wherebeg+10, $where-$wherebeg-10);
      $wherebeg=$where + 1;
    
      # Add colon and backslashes before appending directories
      $newlibpath= $newlibpath . ":\\\\";

      # Find directories and insert in windows' style path
      while ($where >=0) 
      {
        $where = index($ARGV[$i],"/",$wherebeg);
        if ($where < 0) 
        {
          $part = substr($ARGV[$i], $wherebeg, 100);
          $newlibpath = $newlibpath . $part;
        } 
        else 
        {
          $part = substr($ARGV[$i], $wherebeg, $where-$wherebeg);
          $newlibpath = $newlibpath . $part . "\\\\";
        }
        $wherebeg = $where + 1;
      }    
    } 
    else 
    {
      # Assume already have windows path
      $newlibpath = $pathpart;
    }
    $output = $output . "/link /libpath:" . $newlibpath . " ";
  } 
  elsif (index($ARGV[$i],"-g")==0) 
  {
    # Fix the -g problem
    #Do nothing for now -g only generates debugging information for GDB
    # Fix the problem when absolute cygwin paths are used (-I/cygdrive/c/...)
  } 
  elsif (index($ARGV[$i],"/cygdrive")>=0) 
  {
    #Grab the part of the argument before "cygdrive" and preserve it.
    $wherebeg=index($ARGV[$i],"/cygdrive");
    $newpathname=substr($ARGV[$i],0,$wherebeg);

    # Grab the name of the disk, which is expected to be the next directory
    # after cygdrive (ex. /cygdrive/c/Trilinos/..., "c" is the disk name)
    $part = "";
    $where=index($ARGV[$i],"/",$wherebeg+10);
    $newpathname=$newpathname . substr($ARGV[$i], $wherebeg+10, $where-$wherebeg-10);
    $wherebeg=$where + 1;

    # Add colon and backslashes before appending directories
    $newpathname= $newpathname . ":\\\\";

    # Find directories and insert in windows' style path
    while ($where >=0) 
    {
      $where = index($ARGV[$i],"/",$wherebeg);
      if ($where < 0) 
      {
        $part = substr($ARGV[$i], $wherebeg, 100);
        $newpathname = $newpathname . $part;
      } 
      else 
      {
        $part = substr($ARGV[$i], $wherebeg, $where-$wherebeg);
        $newpathname = $newpathname . $part . "\\\\";
      }
      $wherebeg = $where + 1;
    }
    $output = $output . $newpathname . " ";
    #print("$ARGV[$i]","\t","$newpathname","\n");
  }
  elsif (index($ARGV[$i],"-I/home")>=0) 
  {
    #Grab the part of the argument before "home" and preserve it.
    $wherebeg=index($ARGV[$i],"/home");
    $newpathname=substr($ARGV[$i],0,$wherebeg);
          
    # here we have to assume that c:\cygwin is where /home resides
    $newpathname=$newpathname . "c:\\\\cygwin\\\\";
    
    $part = "";
    $where=index($ARGV[$i],"/",$wherebeg);
    $wherebeg=$where + 1;
    
    # Find directories and insert in windows' style path
    while ($where >=0) 
    {
      $where = index($ARGV[$i],"/",$wherebeg);
      if ($where < 0) 
      {
        $part = substr($ARGV[$i], $wherebeg, 100);
        $newpathname = $newpathname . $part;
      } 
      else 
      {
        $part = substr($ARGV[$i], $wherebeg, $where-$wherebeg);
        $newpathname = $newpathname . $part . "\\\\";
      }
      $wherebeg = $where + 1;
    }    
    $output = $output . $newpathname . " ";
    #print("$ARGV[$i]","\t","$newpathname","\n");
    #Otherwise, just pass the argument through unchanged.
  } 
  elsif($ARGV[$i] =~ /\.((c)|(cc)|(cpp)|(c\+\+)|(cxx)|(f)|(f77)|(f90))$/i )
  {
    # found a source file.  Just paste this after the compiler name
    $compiler = $compiler . $ARGV[$i] . " ";
  }
  else 
  {
    $output = $output . $ARGV[$i] . " ";
  }
  $i++;
}

$output = $compiler . $output;
# Let the user know what the command line will be, then execute it.
print ("$output","\n");
exec ("$output") || die "Cannot run reformatted compiler line";

