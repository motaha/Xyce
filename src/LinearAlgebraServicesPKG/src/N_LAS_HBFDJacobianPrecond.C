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
// Filename       : $RCSfile: N_LAS_HBFDJacobianPrecond.C,v $
//
// Purpose        : Implementation file for the Iterative linear solver
//                  interface.
//
// Special Notes  :
//
// Creator        : Heidi K. Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 11/11/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.5 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

// ---------- Standard Includes ----------

#include <Xyce_config.h>

#include <sstream>

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>

#include <N_LAS_HBFDJacobianPrecond.h>
#include <N_LAS_HBFDJacobianEpetraOperator.h>
#include <N_LAS_HBBuilder.h>
#include <N_LOA_HBLoader.h>
#include <N_MPDE_State.h>

#include <N_LAS_MultiVector.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Matrix.h>

#include <N_LAS_Problem.h>
#include <N_LAS_Builder.h>
#include <N_LAS_System.h>
#include <N_LOA_Loader.h>

#include <N_UTL_Timer.h>

#include <N_ERH_ErrorMgr.h>

#include <N_PDS_SerialComm.h>

#include <Teuchos_RCP.hpp>
#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Amesos.h>

using Teuchos::RCP;
using Teuchos::rcp;

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianPrecond::N_LAS_HBFDJacobianPrecond
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
N_LAS_HBFDJacobianPrecond::N_LAS_HBFDJacobianPrecond()
  : N_LAS_Preconditioner()
{
  setDefaultOptions();
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianPrecond::setDefaultOptions
// Purpose       : resets Ifpack options
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBFDJacobianPrecond::setDefaultOptions()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianPrecond::setDefaultOption
// Purpose       : resets Ifpack option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBFDJacobianPrecond::setDefaultOption( const std::string & option )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianPrecond::setOptions
// Purpose       : sets Ifpack options and params from modelblock
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBFDJacobianPrecond::setOptions( const N_UTL_OptionBlock & OB )
{
  // Set the parameters from the list
  std::list<N_UTL_Param>::const_iterator it_tpL = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator end_tpL = OB.getParams().end();
  for (; it_tpL != end_tpL; ++it_tpL)
    {
      this->setParam( *it_tpL );
    }

  // store for restart of solver_
  if( &OB != options_.get() )
    {
      options_ = Teuchos::rcp( &OB, false );
    }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianPrecond::setParam
// Purpose       : sets Ifpack option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBFDJacobianPrecond::setParam( const N_UTL_Param & param )
{
  std::string tag = param.tag();
  std::string uTag = param.uTag();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianPrecond::initGraph
// Purpose       : Initialize the graph using information from the N_LAS_System
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 12/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBFDJacobianPrecond::initGraph( const Teuchos::RCP<N_LAS_Problem> & problem )
{
  // Since this is a time-domain preconditioner, rather than a frequency-domain preconditioner,
  // the graph from the application builder can be used to generate an Amesos solver for each
  // diagonal block.
  RCP<N_LAS_BlockVector> bXt = hbBuilderPtr_->createTimeDomainBlockVector();
  N_ = bXt->blockCount();

  // Generate the vectors for the N_ linear problems to be solved on the diagonal.
  diagRHS_ = rcp( appBuilderPtr_->createVector() );
  diagSoln_ = rcp( appBuilderPtr_->createVector() );
  diagMatrix_.resize(N_);
  epetraProblem_.resize(N_);
  amesosPtr_.resize(N_);

  Amesos amesosFactory;

  for (int i=0; i<N_; ++i) {
    diagMatrix_[i] = rcp( appBuilderPtr_->createMatrix() );
    epetraProblem_[i] = rcp( new Epetra_LinearProblem( &(diagMatrix_[i]->epetraObj()), 
                             &(diagRHS_->epetraObj()), &(diagSoln_->epetraObj()) ) );
    amesosPtr_[i] = rcp( amesosFactory.Create( "Klu", *epetraProblem_[i] ) );
    amesosPtr_[i]->SymbolicFactorization();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianPrecond::initValues
// Purpose       : Initialize the values of the Epetra_LinearSystem's
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 12/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBFDJacobianPrecond::initValues( const Teuchos::RCP<N_LAS_Problem> & problem )
{
  // Get the stored Jacobian matrices from the HB loader.
  std::vector<Teuchos::RCP<N_LAS_Matrix> > vecAppdQdx = hbLoaderPtr_->getStoredQdx();
  std::vector<Teuchos::RCP<N_LAS_Matrix> > vecAppdFdx = hbLoaderPtr_->getStoredFdx();
 
  // Compute the diagonal matrices
  // C(t_i)/h_i + G(t_i)
  double h=0.0;
  for( int i = 0; i < N_; ++i )
  {
    // Initialize the diagonal matrix.
    diagMatrix_[i]->put( 0.0 );

    // If there is only one time step, then h is considered constant.
    if (timeSteps_.size() != 1)
      h = timeSteps_[i];
    else
      h = timeSteps_[0];

    diagMatrix_[i]->linearCombo( 1.0/h, *vecAppdQdx[i], 1.0, *vecAppdFdx[i] );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianPrecond::compute
// Purpose       : Compute a preconditioner M such that M ~= A^{-1}.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_LAS_HBFDJacobianPrecond::compute()
{
  bool precStatus = true;

  // Compute numeric factorization for each block.
  for ( int i=0; i<N_; ++i ) {
    amesosPtr_[i]->NumericFactorization();
  }

  if ( Teuchos::is_null( epetraPrec_ ) )
    epetraPrec_ = fdJacobianOperator( epetraProblem_, amesosPtr_, hbBuilderPtr_, hbLoaderPtr_, timeSteps_ );

  if ( Teuchos::is_null( epetraPrec_ ) )
    return false;

  return precStatus;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBFDJacobianPrecond::apply
// Purpose       : Calls the actual preconditioner to apply y = M*x.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
int N_LAS_HBFDJacobianPrecond::apply( N_LAS_MultiVector & x, N_LAS_MultiVector & y )
{
  int precStatus = 0;
  // If there is no preconditioner to apply return a nonzero code
  if( Teuchos::is_null(epetraPrec_) )
    precStatus = -1;
  else
    precStatus = epetraPrec_->Apply( x.epetraObj(), y.epetraObj() );

  return precStatus;
}

