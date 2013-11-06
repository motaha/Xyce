#!/usr/bin/perl
#
# usage:
#
# cygwin_lib.pl <microsoft lib.exe> <ar args>
#
# This script was originally developed to use lib.exe on cygwin/Windows
# platforms to create libraries.  Unlike the ar command where there is a
# space between the options and the archive name, lib.exe uses the format:
#
# lib.exe /OUT:<archivename> object files...
#
# (1) the first argument passed to this script is the lib.exe program
# to allow for its being installed anywhere
#
# Author : Richard Schiek
# Date : February 21, 2006
#

$MsLib = shift @ARGV;
$cwd = `pwd`;

$i = 0;
$numinputs = @ARGV;

$output = "/usr/bin/lib.exe /OUT:";
$output = "$MsLib";

while ($i < $numinputs) {
  # print "argv[$i] = $ARGV[$i]\n";
  if ( $ARGV[$i]=="x" ) 
  {
    # windows lib.exe doesn't support "ar x" functionality to extract items from
    # an archive.  We need to build a list of archive memebers and then pull
    # them out one at a time.
    $i++;
    @libResponse = `$MsLib /LIST $ARGV[$i]`;
    @libMembers = ();

    # filter out lines that don't end in .obj
    foreach $libMember (@libResponse) 
    { 
      if( $libMember =~ s/\.obj(\r\n)|(\n\r)/\.obj/ )
      {
        push @libMembers, ($libMember);
      }
    }
   
    # now try extracting them 
    foreach $libMember (@libMembers) 
    { 
      `$output /EXTRACT:$libMember $ARGV[$i]; cp $libMember $cwd`;
    }
  }
  elsif( ($ARGV[$i] == "cr") || ($ARGV[$i] == "c") )
  {
    $i++;
    $output = $output . " /OUT:" . $ARGV[$i];
  }
  else
  {
    $filename = $ARGV[$i];
    $filename =~ s/\.a$/\.lib/;    
    $output = $output . $filename . " ";
    $i++;
  }
  $i++;
}

print "Lib command will be: \"$output\"\n";
#exec ("$output") || die "Cannot run reformatted archiving line";


