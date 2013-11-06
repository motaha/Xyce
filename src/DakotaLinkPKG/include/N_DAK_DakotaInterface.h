//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
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
// Filename       : $N_DAK_DakotaInterface.h$
//
// Purpose        : This class defines an interface between the optimization
//                  and analysis routines in Dakota and the circuit routines
//                  in Xyce.  Since this class derives from a Dakota class,
//                  Dakota can work through it to access Xyce.
//
// Special Notes  :
//
//
// Creator        : Richard L. Schiek
//
// Creation Date  : 09/07/2005
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.20.4.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:32 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DAK_DakotaInterface_h
#define Xyce_N_DAK_DakotaInterface_h 1

// Dakota includes
#define DISABLE_DAKOTA_CONFIG_H   1 // Don't need dakota's dakota_config.h 
#include <DirectApplicInterface.H>
#include <ProblemDescDB.H>
#include <DakotaString.H>
#undef DISABLE_DAKOTA_CONFIG_H     // Don't need dakota's dakota_config.h 

// Trilinos includes
#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include <N_UTL_Xyce.h>

// forward declarations
class N_CIR_Xyce;

// standard includes
#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

#include <string>
#include <vector>

class N_DAK_DakotaInterface: public Dakota::DirectApplicInterface
{
public:

#ifdef Xyce_PARALLEL_MPI
    N_DAK_DakotaInterface( const Dakota::ProblemDescDB& problem_db, const MPI_Comm & analysis_comm);
#else
    N_DAK_DakotaInterface( const Dakota::ProblemDescDB& problem_db );
#endif

    ~N_DAK_DakotaInterface();
    
    void setArguments( int iargsIn, char * cargsIn [] );

protected:
    int derived_map_if( const Dakota::String& if_name );
    int derived_map_ac( const Dakota::String& ac_name );
    int derived_map_of( const Dakota::String& of_name );

private:

#ifdef Xyce_PARALLEL_MPI
    // in parallel this is our comm object
    const MPI_Comm MPICommObject;
#endif
    
    // pointers to other objects that N_CIR_Xyce makes
    // but that we'll need to work with
    RCP< N_CIR_Xyce > xyceCirPtr_;
    
    // these are the arguments that Xyce would have received from the command line
    int xyceIargs;
    char ** xyceCargs;

    // This is how we pass variable names and values into Xyce
    vector< pair< string, string > > variableSubVec;
    
    // a Dakota object we'll need to get response variable names
    const Dakota::ProblemDescDB & theProblemDescDB;

    int numResponseVars;

    // The first vector is to hold just then endpoint simulation values
    // when the user has requested that as the response function.  The
    // second vector holds vectors of data for when the user wants to 
    // calculate a norm between a simulated and external signal.
    RCP< std::vector< double > > simulationDataValues;

    void setUpResponse();
    void copyCargs( const int originalIargs, char ** const originalCargs, char ** & copyCargs );
    void deleteCargs( const int len, char ** & cargs );
};


#endif

