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
// Filename       : $RCSfile: N_LAS_TrilinosPrecondFactory.C,v $
//
// Purpose        : Implementation of Trilinos Preconditioning Factory
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 10/01/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.9 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_LAS_TrilinosPrecondFactory.h>
#include <N_LAS_IfpackPrecond.h>
#include <N_LAS_NoPrecond.h>
#ifdef Xyce_ML
#include <N_LAS_MLPrecond.h>
#endif

#include <N_LAS_Problem.h>

#include <N_ERH_ErrorMgr.h>

#include <N_UTL_OptionBlock.h>

// ---------- Trilinos Includes ----------

#include <Epetra_Operator.h>

//-----------------------------------------------------------------------------
// Function      : N_LAS_TrilinosPrecondFactory::N_LAS_TrilinosPrecondFactory
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
N_LAS_TrilinosPrecondFactory::N_LAS_TrilinosPrecondFactory( const N_UTL_OptionBlock & OB )
{
  OB_ = Teuchos::rcp( &OB, false );
  precType_ = "IFPACK";

  std::list<N_UTL_Param>::const_iterator itPI = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator endPI = OB.getParams().end();

  for( ; itPI != endPI; ++itPI )
  {
    if( itPI->uTag() == "PREC_TYPE" && itPI->usVal() != "DEFAULT" )
      precType_ = itPI->usVal();
  }
}
//-----------------------------------------------------------------------------
// Function      : N_LAS_TrilinosPrecondFactory::create
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
Teuchos::RCP<N_LAS_Preconditioner>
N_LAS_TrilinosPrecondFactory::create( const Teuchos::RCP<N_LAS_Problem> & problem )
{
  Teuchos::RCP<N_LAS_Preconditioner> precond;

  // If the problem is matrix free that we can't use any of the available preconditioners.
  if (problem->matrixFree())
  {
    precType_ = "NONE";
  }

  if (precType_ == "IFPACK")
  {
    precond = Teuchos::rcp( new N_LAS_IfpackPrecond() );

    // First set the options for the preconditioner.
    precond->setOptions( *OB_ );

    // Now the preconditioner can create the graph for the preconditioner
    // Note:  The fill and overlap must be set before this point.
    precond->initGraph( problem );
  }
#ifdef Xyce_ML
  else if (precType_ == "ML") {
    precond = Teuchos::rcp( new N_LAS_MLPrecond() );

    // First set the options for the preconditioner.
    precond->setOptions( *OB_ );

    // Now the preconditioner can create the graph for the preconditioner
    // Note:  The fill and overlap must be set before this point.
    precond->initGraph( problem );
  }
#endif
  else if (precType_ == "NONE") {
    // Create an empty preconditioner, which does nothing.
    precond = Teuchos::rcp( new N_LAS_NoPrecond() );
  }
  else
  {
    Xyce::Report::DevelFatal0().in("N_LAS_TrilinosPrecondFactory::create()") << "preconditioning type " << precType_ << " unrecognized!";
  }
  return precond;
}
