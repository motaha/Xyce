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
// Filename       : $RCSfile: N_LAS_SolverFactory.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/12/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.23 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include "Epetra_LinearProblem.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"

// ----------   Xyce Includes   ----------

#include <N_LAS_SolverFactory.h>
#include <N_LAS_Problem.h>

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_SimpleSolver.h>
#include <N_LAS_AmesosSolver.h>
#include <N_LAS_AztecOOSolver.h>
#ifdef Xyce_BELOS
#include <N_LAS_BelosSolver.h>
#endif
#ifdef Xyce_KSPARSE
#include <N_LAS_KSparseSolver.h>
#endif
#ifdef Xyce_SHYLU
#include <N_LAS_ShyLUSolver.h>
#endif

#include <N_IO_CmdParse.h>
#include <N_UTL_OptionBlock.h>

//-----------------------------------------------------------------------------
// Function      : N_LAS_SolverFactory::create
// Purpose       :
// Special Notes : Static
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 05/24/04
//-----------------------------------------------------------------------------
N_LAS_Solver * N_LAS_SolverFactory::create( N_UTL_OptionBlock & options,
                                            N_LAS_Problem & prob,
                                            N_IO_CmdParse & commandLine_ )
{

  int lsDim = (prob.epetraObj().GetRHS())->GlobalLength();

  //If the linear system is trivial, i.e. the matrix is 1x1, then just create a simple solver
  if (lsDim == 1)
    return new N_LAS_SimpleSolver( prob, options );

#ifdef Xyce_PARALLEL_MPI
  std::string type = "AZTECOO";
  if (!prob.matrixFree() && lsDim < 1000)
  {
    type = "KLU";
  }
#else
  std::string type = "KLU";
  if (prob.matrixFree())
  {
    type = "AZTECOO";
  }
#endif

  N_UTL_OptionBlock::ParameterList::iterator itPI = options.begin();
  N_UTL_OptionBlock::ParameterList::iterator endPI = options.end();
  for( ; itPI != endPI; ++itPI )
  {
    if( itPI->uTag() == "TYPE" && itPI->usVal() != "DEFAULT" )
    {
      type = itPI->usVal();
    }
  }

  //Support for resetting linear solver from command line
  ExtendedString CLType = commandLine_.getArgumentValue( "-linsolv" );
  CLType.toUpper();
  if( CLType != "" ) type = CLType;

  // If the linear problem is matrix free, make sure an iterative method is being used.
  if (prob.matrixFree())
  {
    if ((type != "AZTECOO") && (type != "BELOS"))
    {
      std::string msg = "The linear solver option that was specified is not compatible with a matrix free analysis type, changing to AZTECOO";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING, msg);
      type = "AZTECOO";
    }
  }

#ifdef Xyce_PARALLEL_MPI
  // If this is a parallel build and serial execution, use Belos.
  // Otherwise, AztecOO will throw an error because of the serial communicator.
  // Because the N_PDS_ParMap is not guaranteed to be created for any N_LAS_MultiVector, call Epetra.
  int numProcs = ((prob.epetraObj().GetRHS())->Map()).Comm().NumProc();
#ifdef Xyce_BELOS
  if ( numProcs==1 && type=="AZTECOO" )
  {
    type = "BELOS";
  }
#endif
#endif

  if( type == "AZTECOO" )
    return new N_LAS_AztecOOSolver( prob, options );
#ifdef Xyce_BELOS
  else if( type == "BELOS" )
    return new N_LAS_BelosSolver( prob, options );
#endif
#ifdef Xyce_KSPARSE
  else if( type == "KSPARSE" )
    return new N_LAS_KSparseSolver( prob, options );
#endif
#ifdef Xyce_SHYLU
  else if( type == "SHYLU" )
    return new N_LAS_ShyLUSolver( prob, options );
#endif
  else
    return new N_LAS_AmesosSolver( type, prob, options );

  return 0;
}
