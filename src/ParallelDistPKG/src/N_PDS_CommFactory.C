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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_PDS_CommFactory.C,v $
//
// Purpose        : Implementation file for the Comm abstract factory.
//
// Special Notes  : GoF Abstract Factory
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/26/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.13 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_PDS_CommFactory.h>

#include <N_PDS_ParComm.h>
#include <N_PDS_SerialComm.h>
#include <N_PDS_SerialParComm.h>

//-----------------------------------------------------------------------------
// Function      : N_PDS_CommFactory::create
// Purpose       : Creates a serial or parallel comm object based on
//                 system.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_Comm * N_PDS_CommFactory::create( int iargs,
                                        char * cargs[], 
                                        Xyce::Parallel::Machine comm)
{
  N_PDS_Comm * theComm = NULL;
  
#ifdef Xyce_PARALLEL_MPI
  /*
  if( comm )
    return new N_PDS_ParComm( comm );
  else
    return new N_PDS_ParComm( iargs, cargs );
  */
  if( comm )
    theComm = new N_PDS_ParComm( comm );
  else
    theComm = new N_PDS_ParComm( iargs, cargs );
    
  if( theComm->numProc() == 1 )
  {
    // A parallel version of Xyce is trying to run on one processor
    // we'll return a Serial comm object, but store the ParComm
    // object for later deletion.  We do this because other libraries
    // may need MPI running and deleting our N_PDS_ParComm object 
    // would destroy MPI if we created it.
    N_PDS_SerialParComm * serialParComm = new N_PDS_SerialParComm( dynamic_cast<N_PDS_ParComm *>(theComm) );
    theComm = serialParComm;
  }
  return theComm;
#else
  return new N_PDS_SerialComm();
#endif
}

