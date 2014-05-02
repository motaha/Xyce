//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename       : N_UTL_Version.C
//
// Purpose        : set version string  
// Special Notes  : 
//
// Creator        : Eric Rankin
//
// Creation Date  : 
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.19 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>




// ---------- Xyce Includes ---------------------------------------------------
#include <N_UTL_Misc.h>

#include <N_UTL_Version.h>

#include "timestamp.h"

#include <sstream> 

//-----------------------------------------------------------------------------
// Function      : N_UTL_Version::getFullVersionString
// Purpose       : get full banner string for Xyce version
// Special Notes : version should be properly formatted in configure.ac
//               : -Dmacros used to get version and timestamp inserted
// Scope         : 
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
std::string N_UTL_Version::getFullVersionString() 
{
  const std::string tmpVer( VERSION );
  std::string version;
    
  // create developement version string
  if( tmpVer[ 0 ] == 'D' || tmpVer[ 0 ] == 'd' )
  {
    version += "DEVELOPMENT-";
    
    // add the timestamp

    std::ostringstream ver("");

    ver << XYCEBUILDTIMESTAMP;

    version += std::string( ver.str() );
  }

  // create release version string
  else if( tmpVer[ 0 ] == 'R' || tmpVer[ 0 ] == 'r' )
  {
    // prepend release phase if necessary
    if( tmpVer[ 2 ] != ':' )  
    {
      version += "(" + tmpVer.substr( 2, 1 ) + ")";
    }

    // add the release major-minor-rev number 
    int i = tmpVer.find_last_of( ":" );
    version += "Release " + tmpVer.substr( i + 1, tmpVer.length() - i );
  }

  version += getBuildVariant();

  return version;
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_Version::getShortVersionString
// Purpose       : get the maj-min-rev number for Xyce version
// Special Notes : version should be properly formatted in configure.ac
// Scope         : 
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
std::string N_UTL_Version::getShortVersionString()
{
  const std::string tmpVer( VERSION );

  // get position of the major-minor-rev number  
  int i = tmpVer.find_last_of( ":" );
  
  return tmpVer.substr( i + 1, tmpVer.length() - i );
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_Version::getBuildVariant
// Purpose       : Return a short string indicating special builds
// Special Notes : Current special builds are:
//               :  (blank) = full ECI build
//               :  OS      = open source build
//               :  norad   = no rad models, but not open source
// Scope         : 
// Creator       : Tom Russo
// Creation Date : 6/10/2013
//-----------------------------------------------------------------------------
std::string N_UTL_Version::getBuildVariant()
{
  std::string variant;
#ifdef Xyce_RAD_MODELS
  variant="";
#else
#ifdef Xyce_NONFREE_MODELS
  variant="-norad";
#else
  variant="-opensource";
#endif // Xyce_NONFREE_MODELS
#endif // Xyce_RAD_MODELS

  return variant;
}


  
//-----------------------------------------------------------------------------
// Function      : N_UTL_Version::getCapabilities
// Purpose       : Return a string indicating compiled-in features
// Special Notes : 
// Scope         : 
// Creator       : Tom Russo
// Creation Date : 6/10/2013
//-----------------------------------------------------------------------------
std::string N_UTL_Version::getCapabilities()
{
  std::string capabilities="";

#ifdef Xyce_PARALLEL_MPI
  capabilities += "Parallel with MPI\n";
#else
  capabilities += "Serial\n";
#endif

#ifdef Xyce_RAD_MODELS
  capabilities += "Radiation models\n";
#endif

#ifdef Xyce_NONFREE_MODELS
  capabilities += "Non-GPL device models\n";
#endif

#ifdef Xyce_USE_FFT
  capabilities += "FFT";
#ifdef Xyce_USE_INTEL_FFT
  capabilities += "(Intel FFT)\n";
#else
  capabilities += "(FFTW)\n";
#endif
#endif

#ifdef Xyce_USE_HDF5
  capabilities += "HDF5\n";
#endif

#ifdef Xyce_REACTION_PARSER
  capabilities += "Reaction parser\n";
#endif

#ifdef Xyce_ATHENA
  capabilities += "ATHENA\n";
#endif

#ifdef Xyce_VERBOSE_TIME
  capabilities += "Verbose output - time integrator\n";
#endif

#ifdef Xyce_VERBOSE_LINEAR
  capabilities += "Verbose output - linear solver\n";
#endif

#ifdef Xyce_VERBOSE_NONLINEAR
  capabilities += "Verbose output - nonlinear solver\n";
#endif

  return capabilities;
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_Version::getLicense
// Purpose       : Return a string indicating license
// Special Notes : 
// Scope         : 
// Creator       : Tom Russo
// Creation Date : 6/10/2013
//-----------------------------------------------------------------------------
std::string N_UTL_Version::getLicense()
{
  std::string License="";

#ifdef Xyce_RAD_MODELS
  License+= " \nEXPORT CONTROLLED SOFTWARE\n\n";
  License+="Copyright (c) 2003, Sandia Corporation, Albuquerque, NM, USA.  Under the\n";
  License+= " terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for\n";
  License+= " use of this work by or on behalf of the U.S. Government.  Export of this\n";
  License += " program may require a license from the United States Government.\n";
#else
#ifdef Xyce_NONFREE_MODELS
  License = "THIS IS NOT AN OPEN SOURCE BINARY.\n";

  License +="The EKV3 model is not an open source model.  It was developed by the\n";
  License +="EKV Team of the Electronics Laboratory-TUC of the Technical University\n";
  License +="of Crete.  TUC has granted Sandia National Laboratories a\n";
  License +="non-exclusive, roalty-free licnse to use, copy, modify, implement and\n";
  License +="distribute the EKV3 Model Code in the form of compiled, executable\n";
  License +="code.  Sandia National Laboratories is not licensed to distribute this\n";
  License +="code in source form.  Documentation of the EKV3 model is available on\n";
  License +="the Xyce web site, http://charleston.sandia.gov/xyce/xyce_internal.html, and more information about the\n";
  License +="EKV model (including contact information for the EKV Model Team) is\n";
  License +="available at the EKV web site, http://ekv.epfl.ch/.\n";
  License +="\n";
  License +="All components of this version of Xyce OTHER THAN the EKV3 model are\n";
  License +="released under the GNU Public License:\n";
#endif

  License +="\n    Xyce(TM) Parallel Electrical Simulator\n\n";
  License +="    Copyright 2002-2014 Sandia Corporation. Under the terms\n";
  License +="    of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.\n";
  License +="    Government retains certain rights in this software.\n";
  License +="\n";
  License +="    This program is free software: you can redistribute it and/or modify\n";
  License +="    it under the terms of the GNU General Public License as published by\n";
  License +="    the Free Software Foundation, either version 3 of the License, or\n";
  License +="    (at your option) any later version.\n";
  License +="\n";
  License +="    This program is distributed in the hope that it will be useful,\n";
  License +="    but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
  License +="    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
  License +="    GNU General Public License for more details.\n";
  License +="\n";
  License +="    You should have received a copy of the GNU General Public License\n";
  License +="    along with this program.  If not, see <http://www.gnu.org/licenses/>.\n";
  License +="\n";
  License +="\n";

#ifdef Xyce_NONFREE_MODELS
  License +="To obtain source code for components of Xyce EXCLUDING the EKV model,\n";
  License +="see the Xyce web site, http://xyce.sandia.gov/.\n";
#endif

#endif

  return License;
}
