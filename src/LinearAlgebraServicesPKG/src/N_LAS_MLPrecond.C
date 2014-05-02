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
// Filename       : $RCSfile: N_LAS_MLPrecond.C,v $
//
// Purpose        : Implementation file for the Iterative linear solver
//                  interface.
//
// Special Notes  :
//
// Creator        : Heidi K. Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 09/27/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#include <sstream>

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>

#include <N_LAS_Problem.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_MLPrecond.h>

#include <N_UTL_OptionBlock.h>

#include <N_UTL_Timer.h>

#include <N_ERH_ErrorMgr.h>

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <ml_MultiLevelPreconditioner.h>



//-----------------------------------------------------------------------------
// Function      : N_LAS_MLPrecond::N_LAS_MLPrecond
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
N_LAS_MLPrecond::N_LAS_MLPrecond()
  : N_LAS_Preconditioner()
{
  setDefaultOptions();
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MLPrecond::setDefaultOptions
// Purpose       : resets ML options
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_MLPrecond::setDefaultOptions()
{
  ML_Epetra::SetDefaults("SA",mlList_);
  maxLevel_ = maxLevel_default_;
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MLPrecond::setDefaultOption
// Purpose       : resets ML option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_MLPrecond::setDefaultOption( const string & option )
{
  if (option == "ML_MAX_LEVEL") {
    maxLevel_ = maxLevel_default_;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MLPrecond::setOptions
// Purpose       : sets ML options and params from modelblock
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_MLPrecond::setOptions( const N_UTL_OptionBlock & OB )
{
  for( list<N_UTL_Param>::const_iterator it_tpL = OB.params.begin();
         it_tpL != OB.params.end(); ++it_tpL )
  {
    string tag = it_tpL->uTag();

    if( tag == "ML_MAX_LEVEL" ) maxLevel_ = it_tpL->iVal();
  }

  // store for restart of solver_
  if( &OB != options_.get() )
    {
      options_ = Teuchos::rcp( new N_UTL_OptionBlock(OB) );
    }
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MLPrecond::setParam
// Purpose       : sets ML option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_MLPrecond::setParam( const N_UTL_Param & param )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MLPrecond::initGraph
// Purpose       : Set up the graph pattern for the preconditioner.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_MLPrecond::initGraph( const Teuchos::RCP<N_LAS_Problem> & problem )
{
  // output level, 0 being silent and 10 verbose
  mlList_.set("output", 0);
  // maximum number of levels
  mlList_.set("max levels", maxLevel_);
  // set finest level to 0
  mlList_.set("increasing or decreasing","increasing");

  // use Uncoupled scheme to create the aggregate
  mlList_.set("aggregation: type", "Uncoupled");

  // use Gauss-Seidel smoother
  mlList_.set("smoother: type","Gauss-Seidel");

  // use both pre and post smoothing
  mlList_.set("smoother: pre or post", "both");

  // solve with serial direct solver KLU
  mlList_.set("coarse: type","Amesos-KLU");

  // Create the precondtioning object, but don't compute it.
  Epetra_CrsMatrix * epetraA = dynamic_cast<Epetra_CrsMatrix*>(problem->epetraObj().GetMatrix());
  mlPrecond_ = Teuchos::rcp( new ML_Epetra::MultiLevelPreconditioner(*epetraA, mlList_, false) );  
  epetraPrec_ = mlPrecond_;  

  return true;  
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MLPrecond::initValues
// Purpose       : Set the values for the preconditioner.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_MLPrecond::initValues( const Teuchos::RCP<N_LAS_Problem> & problem )
{
  bool precStatus = true;
  problem_ = problem;

  return precStatus;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MLPrecond::compute
// Purpose       : Compute a preconditioner M such that M ~= A^{-1}.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_MLPrecond::compute()
{
  bool precStatus = true;
  if ( Teuchos::is_null( epetraPrec_ ) )
    return false;

  // Now compute the preconditioner.
  mlPrecond_->ComputePreconditioner();

  return precStatus;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_MLPrecond::apply
// Purpose       : Calls the actual preconditioner to apply y = M*x.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
int N_LAS_MLPrecond::apply( N_LAS_MultiVector & x, N_LAS_MultiVector & y )
{
  int precStatus = 0;

  // If there is no preconditioner to apply return a nonzero code
  if( Teuchos::is_null(epetraPrec_) )
    precStatus = -1;
  else
    precStatus = epetraPrec_->Apply( x.epetraObj(), y.epetraObj() );

  return precStatus;
}

