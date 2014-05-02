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
// Filename       : Xyce.C
//
// Purpose        : front end for standalone Xyce executable
//
// Special Notes  :
//
// Creator        : Eric Rankin
//
// Creation Date  : 01/28/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.18 $
//
// Revision Date  : $Date: 2014/02/24 23:49:13 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#include <string>

#ifdef HAVE_LINUX_EXCEPTIONS
#include <fenv.h>
#endif

#include <N_CIR_Xyce.h>
#include <N_ERH_ErrorMgr.h>

#ifdef Xyce_Dakota
#include <N_DAK_DakotaController.h>
#endif

// Function to be called if memory runs out:
void _new_handler (void)
{
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, "OUT OF MEMORY (error in 'new')");
}


//
//-----------------------------------------------------------------------------
// Function      : checkForDakotaFlag()
// Purpose       : This function scans the argument line to determine if this
//                 is a Dakota controlled run. (i.e. has -dakota on arg list)
// Special Notes :
// Scope         :
// Creator       : Rich Schiek
// Creation Date : 10/14/2008
//-----------------------------------------------------------------------------
bool checkForDakotaFlag( const int iargs, const char * const cargs[])
{
  bool result = false;

  for( int i=0; i<iargs; i++ )
  {
    std::string currentArg(cargs[i]);
    if( currentArg == "-dakota" )
    {
      result = true;
      break;
    }
  }
  return result;
}


//
//-----------------------------------------------------------------------------
// Function      : main
// Purpose       : front end for standalone Xyce executable
// Special Notes :
// Scope         :
// Creator       : Eric Rankin
// Creation Date : 01/28/2004
//-----------------------------------------------------------------------------
int main( int iargs, char *cargs[] )
{
  // Set divide by zero, and invalid operation handling on linux
#ifdef HAVE_LINUX_EXCEPTIONS
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
  // Set out of memory detection on all systems
  std::set_new_handler (&_new_handler);

  bool bsuccess = false;
#ifdef Xyce_Dakota
  if( checkForDakotaFlag( iargs, cargs ) == true )
  {
    // this will be a Xyce simulation where Dakota creates and
    // then runs N_CIR_Xyce(). So pass control to it

    N_DAK_DakotaController dakotaController( iargs, cargs );
    bsuccess = dakotaController.run();
  }
  else
#endif
  {
    N_CIR_Xyce xyce;
    bsuccess = xyce.run( iargs, cargs );
  }
  (bsuccess) ? exit(0) : exit(-1);
}

