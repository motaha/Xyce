#!/usr/bin/perl
#
# usage:
#
# mac_osx_ar.pl <ar args>
#
# This script was developed to use "ar" on Mac OS X through libtool.
#
# Under Mac OS X, programs link dynamically to shared libraries. 
# If both a shared and static library are available, the linking
# will be dynamic even if one used flags like "-static".
# 
# A way around this defect is to specify specific library files
# such as /usr/lib/libm.a rather than -lm.  Thus, the intel
# compilers will specify a set of libraries to pass down to the
# linker as /opt/intel/lib/libblah.a
# 
# When the build tools probe a compiler to check for library
# dependencies, these fully specified library names get captured
# and then libtool tries to add them to its own archives.  Normally
# this wouldn't be a problem, but under Mac OS X, some libraries
# can max PPC and Intel x86 code (i.e. universal libraries) and
# when libtool tries to use "ar" to make a normal library  with
# some universal components, this code mixing generates 
# a failure in "ar".
# 
# This script tries to work around these issues by acting as a
# conduit between libtool and "ar"  It scans the arguments passed
# to it and removes any archive files (lib*.a) that do not reside in
# the .libs/ directory that libtool makes.
#
# Author : Richard Schiek
# Date : October 4, 2007
#

$i = 0;
$numinputs = @ARGV;

$output = "ar ";

while ($i < $numinputs) {
  if ( $ARGV[$i] =~ /.*\/\.libs\/.*a$/ )
  {
    $output = $output . $ARGV[$i] . " ";
  }
  elsif ( $ARGV[$i] =~ /^\/.*\.a$/ ) 
  {
    print "Discarding $ARGV[$i]\n";
  }
  else
  {  
    $output = $output . $ARGV[$i] . " ";
  }
  $i++;
}

print "Lib command will be: \"$output\"\n";
exec ("$output") || die "Cannot run reformatted archiving line";



