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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_LAS_HBPrecondFactory.C,v $
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
// Revision Number: $Revision: 1.1.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:45 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_LAS_Preconditioner.h>
#include <N_LAS_HBPrecondFactory.h>
#include <N_LAS_HBBlockJacobiPrecond.h>

#include <N_LAS_Problem.h>
#include <N_LAS_System.h>
#include <N_LAS_NoPrecond.h>

#include <N_ERH_ErrorMgr.h>

// ---------- Trilinos Includes ----------

#include <Epetra_Operator.h>

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBPrecondFactory::N_LAS_HBPrecondFactory
// Purpose       :
// Special Notes : 
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
N_LAS_HBPrecondFactory::N_LAS_HBPrecondFactory()
{
  precType_ = "NONE";
}
//-----------------------------------------------------------------------------
// Function      : N_LAS_HBPrecondFactory::N_LAS_HBPrecondFactory
// Purpose       :
// Special Notes : 
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
N_LAS_HBPrecondFactory::N_LAS_HBPrecondFactory( const N_UTL_OptionBlock & OB )
{
  OB_ = Teuchos::rcp( new N_UTL_OptionBlock( OB ) );
  precType_ = "NONE";

  list<N_UTL_Param>::const_iterator it_tpL = OB.getParams().begin();
  list<N_UTL_Param>::const_iterator end_tpL = OB.getParams().end();
  for( ; it_tpL != end_tpL; ++it_tpL )
  {
    if( it_tpL->uTag() == "PREC_TYPE" && it_tpL->usVal() != "DEFAULT" )
      precType_ = it_tpL->usVal();
  }
}
//-----------------------------------------------------------------------------
// Function      : N_LAS_HBPrecondFactory::create
// Purpose       :
// Special Notes : 
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 10/01/07
//-----------------------------------------------------------------------------
Teuchos::RCP<N_LAS_Preconditioner> 
N_LAS_HBPrecondFactory::create( const Teuchos::RCP<N_LAS_System> & lasSysPtr )
{
  lasSysPtr_ = lasSysPtr;

  Teuchos::RCP<N_LAS_Preconditioner> precond;

  if (precType_ == "NONE") {
    // Create an empty preconditioner, which does nothing.
    precond = Teuchos::rcp( new N_LAS_NoPrecond() );
  }
  else if (precType_ == "BLOCK_JACOBI") {
    precond = Teuchos::rcp( new N_LAS_HBBlockJacobiPrecond() );
    precond->setOptions( *OB_ );
    
    // Register necessary classes for block Jacobi preconditioner.
    Teuchos::RCP<N_LAS_HBBlockJacobiPrecond> tmpPrecond 
      = Teuchos::rcp_dynamic_cast<N_LAS_HBBlockJacobiPrecond>( precond );

    tmpPrecond->registerLinearSystem( lasSysPtr_ );
    tmpPrecond->registerAppLoader( appLoaderPtr_ );
    tmpPrecond->registerAppBuilder( appBuilderPtr_ );
    tmpPrecond->registerHBLoader( hbLoaderPtr_ );
    tmpPrecond->registerHBBuilder( hbBuilderPtr_ );
    tmpPrecond->registerMPDEState( statePtr_ );
    tmpPrecond->registerDeviceInterface( devInterfacePtr_ );
    tmpPrecond->setFastTimes( times_ );

    // Initialize the graph for the preconditioner
    Teuchos::RCP<N_LAS_Problem> tmpProblem;
    tmpPrecond->initGraph( tmpProblem );
  }
  else {
    string msg = "N_LAS_HBPrecondFactory::create()";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR_0, msg + ", preconditioning type " + precType_ + " unrecognized!\n");
  }

  return precond;
}
