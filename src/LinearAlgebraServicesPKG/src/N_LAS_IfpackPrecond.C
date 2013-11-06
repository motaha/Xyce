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
// Filename       : $RCSfile: N_LAS_IfpackPrecond.C,v $
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
// Revision Number: $Revision: 1.19.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:45 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#include <sstream>

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>

#include <N_LAS_IfpackPrecond.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Problem.h>

#include <N_UTL_OptionBlock.h>

#include <N_UTL_Timer.h>

#include <N_ERH_ErrorMgr.h>

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <Ifpack_IlukGraph.h>
#include <Ifpack_CrsRiluk.h>
#include <Ifpack.h>
#include <Ifpack_Preconditioner.h>

// static class member initializations
// Default preconditioner values
const bool   N_LAS_IfpackPrecond::useFactory_default_ = false;
const string N_LAS_IfpackPrecond::ifpackType_default_ = "Amesos";
const double N_LAS_IfpackPrecond::diagPerturb_default_= 0.0;
const int    N_LAS_IfpackPrecond::overlap_default_    = 0;
const double N_LAS_IfpackPrecond::dropTol_default_    = 1.0e-03;
const double N_LAS_IfpackPrecond::ilutFill_default_   = 2.0;
const double N_LAS_IfpackPrecond::rThresh_default_    = 1.0001;
const double N_LAS_IfpackPrecond::aThresh_default_    = 0.0001;

//-----------------------------------------------------------------------------
// Function      : N_LAS_IfpackPrecond::N_LAS_IfpackPrecond
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
N_LAS_IfpackPrecond::N_LAS_IfpackPrecond()
  : N_LAS_Preconditioner()
{
  setDefaultOptions();
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_IfpackPrecond::setDefaultOptions
// Purpose       : resets Ifpack options
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_IfpackPrecond::setDefaultOptions()
{
  // Set defaults
  useFactory_ = useFactory_default_;
  ifpackType_ = ifpackType_default_;
  diagPerturb_ = diagPerturb_default_;
  dropTol_ = dropTol_default_;
  ilutFill_ = ilutFill_default_;
  rThresh_ = rThresh_default_;
  aThresh_ = aThresh_default_;
  overlap_ = overlap_default_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_IfpackPrecond::setDefaultOption
// Purpose       : resets Ifpack option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_IfpackPrecond::setDefaultOption( const string & option )
{
  if( option == "AZ_athresh" )         aThresh_ = aThresh_default_;
  if( option == "AZ_rthresh" )         rThresh_ = rThresh_default_;
  if( option == "AZ_ilut_fill" )       ilutFill_ = ilutFill_default_;
  if( option == "AZ_drop" )            dropTol_ = dropTol_default_;
  if( option == "AZ_overlap" )         overlap_ = overlap_default_;
  if( option == "use_ifpack_factory" ) useFactory_ = useFactory_default_;
  if( option == "ifpack_type" )        ifpackType_ = ifpackType_default_;
  if( option == "diag_perturb" )       diagPerturb_ = diagPerturb_default_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_IfpackPrecond::setOptions
// Purpose       : sets Ifpack options and params from modelblock
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_IfpackPrecond::setOptions( const N_UTL_OptionBlock & OB )
{
  // Set the parameters from the list
  list<N_UTL_Param>::const_iterator it_tpL = OB.getParams().begin();
  list<N_UTL_Param>::const_iterator end_tpL = OB.getParams().end();
  for (; it_tpL != end_tpL; ++it_tpL)
    {
      setParam( *it_tpL );
    }

  // store for restart of solver_
  if( &OB != options_.get() )
    {
      options_ = Teuchos::rcp( new N_UTL_OptionBlock(OB) );
    }

  return STATUS_SUCCESS;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_IfpackPrecond::setParam
// Purpose       : sets Ifpack option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_IfpackPrecond::setParam( const N_UTL_Param & param )
{
  string tag = param.tag();
  string uTag = param.uTag();

  // Set our copies of these parameters that get passed to the solver in the
  // "iterate" command
  if( tag == "AZ_overlap" )
    overlap_ = param.iVal();
  else if( tag == "AZ_athresh")
    aThresh_ = param.dVal();
  else if( tag == "AZ_rthresh")
    rThresh_ = param.dVal();
  else if( tag == "AZ_drop")
    dropTol_ = param.dVal();
  else if( tag == "AZ_ilut_fill")
    ilutFill_ = param.dVal();
  else if( tag == "use_ifpack_factory")
    useFactory_ = param.iVal();
  else if( tag == "ifpack_type")
    ifpackType_ = param.usVal();
  else if( tag == "diag_perturb")
    diagPerturb_ = param.dVal();
  else
    return false;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_IfpackPrecond::initGraph
// Purpose       : Set up the graph pattern for the preconditioner.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_IfpackPrecond::initGraph( const Teuchos::RCP<N_LAS_Problem> & problem )
{
  bool precStatus = true;

  if (useFactory_)
  {
    // Because a new graph is being used, recreate the preconditioner.

    // Create Ifpack factory.
    Ifpack Factory;

    // Create the preconditioner.
    Epetra_CrsMatrix * epetraA = dynamic_cast<Epetra_CrsMatrix*>(problem->epetraObj().GetMatrix());
    ifpackPrecond_ = Teuchos::rcp( Factory.Create(ifpackType_, epetraA, overlap_) );

    if (ifpackPrecond_ == Teuchos::null) {
      std::string msg = "N_LAS_IfpackPrecond::initGraph()";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_ERROR_0, msg + ", preconditioning type " + ifpackType_ + " unrecognized!\n");
    }

    // Pass solver type to Amesos.
    Teuchos::ParameterList List;
    if (ifpackType_ == "Amesos") {
      List.set("amesos: solver type", "Amesos_Klu");
      if (diagPerturb_ != 0.0)
        List.set("AddToDiag", diagPerturb_);
    }
    else {
      List.set("fact: absolute threshold", aThresh_);
      List.set("fact: relative threshold", rThresh_);
      List.set("fact: level-of-fill", (int)ilutFill_);
      List.set("fact: drop tolerance", dropTol_);
    }

    // Set the parameters
    IFPACK_CHK_ERR(ifpackPrecond_->SetParameters(List));

    // Compute symbolic factorization stage of preconditioner.
    IFPACK_CHK_ERR(ifpackPrecond_->Initialize());

    // Set the Epetra object to point to this preconditioner.
    epetraPrec_ = ifpackPrecond_;
  }
  else
  {
    // Create the graph.
    Epetra_CrsMatrix * epetraA = dynamic_cast<Epetra_CrsMatrix*>(problem->epetraObj().GetMatrix());
    const Epetra_CrsGraph & Graph = epetraA->Graph();
    ilukGraph_ = Teuchos::rcp( new Ifpack_IlukGraph( Graph, static_cast<int>(ilutFill_), overlap_ ) );
    int graphRet = ilukGraph_->ConstructFilledGraph();
    if ( graphRet != 0 )
      precStatus = false;

    // Because a new graph was created, destroy the preconditioner to make sure it's recreated.
    rILUK_ = Teuchos::null;
    epetraPrec_ = Teuchos::null;
  }

  return precStatus;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_IfpackPrecond::initValues
// Purpose       : Set the values for the preconditioner.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_IfpackPrecond::initValues( const Teuchos::RCP<N_LAS_Problem> & problem )
{
  bool precStatus = true;
  problem_ = problem;

  if (useFactory_)
  {
    // It is currently assumed that the same N_LAS_Matrix object you are passing in
    // is the one you used to set the graph.  This method takes a pointer to
    // it and will use it during the preconditioner computation.
    if ( Teuchos::is_null( ifpackPrecond_ ) )
      initGraph( problem_ );
  }
  else
  {
    // Initialize the graph if we need to.
    if ( Teuchos::is_null( ilukGraph_ ) )
      initGraph( problem_ );

    // Create the preconditioner if one doesn't exist.
    if ( Teuchos::is_null( rILUK_ ) )
    {
      rILUK_ = Teuchos::rcp( new Ifpack_CrsRiluk( *ilukGraph_ ) );
      rILUK_->SetAbsoluteThreshold( aThresh_ );
      rILUK_->SetRelativeThreshold( rThresh_ );
    }

    Epetra_CrsMatrix * epetraA = dynamic_cast<Epetra_CrsMatrix*>(problem_->epetraObj().GetMatrix());
    int initErr = rILUK_->InitValues( *epetraA );
    if (initErr < 0)
      precStatus = false;

    epetraPrec_ = rILUK_;
  }

  return precStatus;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_IfpackPrecond::compute
// Purpose       : Compute a preconditioner M such that M ~= A^{-1}.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
bool N_LAS_IfpackPrecond::compute()
{
  bool precStatus = true;
  if ( Teuchos::is_null( epetraPrec_ ) )
    return false;

  if (useFactory_)
  {
    // Build the preconditioner using the values in problem_; numeric factorization.
    IFPACK_CHK_ERR(ifpackPrecond_->Compute());
  }
  else
  {
    Epetra_CrsMatrix * epetraA = dynamic_cast<Epetra_CrsMatrix*>(problem_->epetraObj().GetMatrix());
    bool transpose = epetraA->UseTranspose();

    int factErr = rILUK_->Factor();
    if (factErr < 0)
      return false;

#ifdef Xyce_VERBOSE_LINEAR
    double condest;
    rILUK_->Condest( transpose, condest );
#endif

    // Define label for printing out during the solve phase
    ostringstream ost;
    ost << "Ifpack_CrsRiluk Preconditioner: LevelFill = " << ilutFill_ << endl <<
           "                                Overlap = " << overlap_ << endl <<
           "                                Athresh = " << aThresh_ << endl <<
           "                                Rthresh = " << rThresh_ << endl <<
#ifdef Xyce_VERBOSE_LINEAR
           "                                CondEst = " << condest  << endl <<
#endif
           "                                ErrCode = " << factErr  << endl;
    string label = ost.str();
    rILUK_->SetLabel(label.c_str());

    rILUK_->SetUseTranspose( transpose );
  }

  return precStatus;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_IfpackPrecond::apply
// Purpose       : Calls the actual preconditioner to apply y = M*x.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/27/07
//-----------------------------------------------------------------------------
int N_LAS_IfpackPrecond::apply( N_LAS_MultiVector & x, N_LAS_MultiVector & y )
{
  int precStatus = 0;

  // If there is no preconditioner to apply return a nonzero code
  if( Teuchos::is_null(epetraPrec_) )
    precStatus = -1;
  else
    precStatus = epetraPrec_->ApplyInverse( x.epetraObj(), y.epetraObj() );

  return precStatus;
}

