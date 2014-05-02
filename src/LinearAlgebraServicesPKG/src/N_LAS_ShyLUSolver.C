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
// Filename       : $RCSfile: N_LAS_ShyLUSolver.C,v $
//
// Purpose        : Implementation file for the ShyLU linear solver interface.
//
// Special Notes  :
//
// Creator        : Heidi K. Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 09/25/2007
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.20 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>


// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>

#include <N_LAS_ShyLUSolver.h>

#include <N_UTL_OptionBlock.h>

#include <N_UTL_Timer.h>

#include <N_ERH_ErrorMgr.h>

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>

#include <N_LAS_TransformTool.h>
#include <N_LAS_Problem.h>

#include <N_LAS_Matrix.h>

#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_BlockMapOut.h>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Utils.hpp>

//-----------------------------------------------------------------------------
// Function      : N_LAS_ShyLUSolver::~N_LAS_ShyLUSolver
// Purpose       : Default destructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
N_LAS_ShyLUSolver::~N_LAS_ShyLUSolver()
{
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_ShyLUSolver::N_LAS_ShyLUSolver
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
N_LAS_ShyLUSolver::N_LAS_ShyLUSolver( N_LAS_Problem & problem,
				      N_UTL_OptionBlock & options )
: N_LAS_Solver(true),
  symmetryDefault_(1),
  innerMaxIterDefault_(30),
  innerTolDefault_(1.0e-12),
  relThreshDefault_(1e-3),
  diagFactorDefault_(0.05),
  outerSolverDefault_("Belos"),
  separatorTypeDefault_("Wide"),
  schurSolverDefault_("AztecOO-Exact"),  // "AztecOO-Inexact", "Amesos", "Guided Probing"
  schurApproxTypeDefault_("Threshold"),
  outputLS_(0),
  outputBaseLS_(0),
  lasProblem_(problem),
  problem_(&(problem.epetraObj())),
  updatedParams_(false),
  linearResidual_(1.0),
  tProblem_(0)
{
  options_ = Teuchos::rcp( new N_UTL_OptionBlock( options ) );
  timer_ = Teuchos::rcp( new N_UTL_Timer( problem_->GetLHS()->Comm() ) );

  shyluParams_ = Teuchos::rcp( new Teuchos::ParameterList() );
  setDefaultOptions();
  setOptions( *options_ );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_ShyLUSolver::setDefaultOptions
// Purpose       : resets Aztec options
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 9/25/07
//-----------------------------------------------------------------------------
bool N_LAS_ShyLUSolver::setDefaultOptions()
{
  symmetry_ = symmetryDefault_;
  innerMaxIter_ = innerMaxIterDefault_;
  innerTol_ = innerTolDefault_;
  relThresh_ = relThreshDefault_;
  diagFactor_ = diagFactorDefault_;
  outerSolver_ = outerSolverDefault_;
  separatorType_ = separatorTypeDefault_;
  schurSolver_ = schurSolverDefault_;
  schurApproxType_ = schurApproxTypeDefault_;
  updatedParams_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_ShyLUSolver::setDefaultOption
// Purpose       : resets Aztec option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 9/25/07
//-----------------------------------------------------------------------------
bool N_LAS_ShyLUSolver::setDefaultOption( const string & option )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_ShyLUSolver::setOptions
// Purpose       : sets Aztec options and params from modelblock
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool N_LAS_ShyLUSolver::setOptions(const N_UTL_OptionBlock& OB)
{
  shyluParams_->set("Symmetry", symmetry_); 
  shyluParams_->set("Inner Solver MaxIters", innerMaxIter_); 
  shyluParams_->set("Inner Solver Tolerance", innerTol_);
  shyluParams_->set("Relative Threshold", relThresh_);
  shyluParams_->set("Diagonal Factor", diagFactor_); 
  shyluParams_->set("Outer Solver Library", outerSolver_); 
  shyluParams_->set("Separator Type", separatorType_); 
  shyluParams_->set("Schur Approximation Method", schurApproxType_); 
  shyluParams_->set("Schur Complement Solver", schurSolver_); 
  //shyluParams_->set("Schur Recompute Iteration", 5);

  list<N_UTL_Param>::const_iterator it_tpL = OB.getParams().begin();
  list<N_UTL_Param>::const_iterator end_tpL = OB.getParams().end();
  for (; it_tpL != end_tpL; ++it_tpL)
  {
    setParam( *it_tpL );
  }
  
  // store for restart of solver_
  if( &OB != &*options_ )
  {
    options_ = Teuchos::rcp( new N_UTL_OptionBlock(OB) );
  }

  if (!lasProblem_.matrixFree()) 
  {
    // create the transformation object if needed.
    if( Teuchos::is_null(transform_) ) transform_ = N_LAS_TransformTool()( OB );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_ShyLUSolver::setParam
// Purpose       : sets Aztec option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool N_LAS_ShyLUSolver::setParam( const N_UTL_Param & param )
{
  string tag = param.tag();
  string uTag = param.uTag();

  // Set our copies of these parameters.
  if( tag == "AZ_max_iter" )
    setInnerMaxIter(param.getImmutableValue<int>());
  else if( tag == "AZ_tol" )
    setInnerTol(param.getImmutableValue<double>());
  else if( tag == "ShyLU_rthresh" )
    setRelThresh(param.getImmutableValue<double>());
  else if( uTag == "OUTPUT_LS" )
    outputLS_ = param.getImmutableValue<int>();
  else if( uTag == "OUTPUT_BASE_LS" )
    outputBaseLS_ = param.getImmutableValue<int>();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_ShyLUSolver::getInfo
// Purpose       : sets Aztec option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool N_LAS_ShyLUSolver::getInfo( N_UTL_Param & param )
{
  if( param.tag() == "AZ_max_iter" )
    param.setVal( innerMaxIter_ );
  else if( param.tag() == "Iterations" )
    param.setVal( (int)numLinearIters_ );
  else if( param.tag() == "AZ_tol" )
    param.setVal( innerTol_ );
  else 
    return false;
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_ShyLUSolver::setShyLUOption
// Purpose       : sets ShyLU option
// Special Notes : Takes a string as the option identifier
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Computational Sciences
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool N_LAS_ShyLUSolver::setShyLUOption_(const char * paramName,
                                            const int val)
{
  return setShyLUCntl_( N_UTL_Param(paramName, val) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_ShyLUSolver::setShyLUParam
// Purpose       : sets ShyLU parameter
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool N_LAS_ShyLUSolver::setShyLUParam_(const char * paramName,
                                         const double val)
{
  return setShyLUCntl_( N_UTL_Param(paramName, val) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_ShyLU::printParams_
// Purpose       : Print out the linear solver parameter values.
// Special Notes :
// Scope         : Private
// Creator       : Heidi Thornquist, SNL, Computational Sciences
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
void N_LAS_ShyLUSolver::printParams_() const
{}

//-----------------------------------------------------------------------------
// Function      : N_LAS_ShyLUSolver::solve
// Purpose       : Calls the actual solver to solve Ax=b.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
int N_LAS_ShyLUSolver::solve( bool ReuseFactors )
{
  // Start the timer...
  timer_->resetStartTime();

#if defined(Xyce_DETAILED_LINSOLVE_TIMES) || defined(Xyce_VERBOSE_LINEAR)
  double time1 = 0.0, time2 = 0.0;
#endif

#ifdef Xyce_DETAILED_LINSOLVE_TIMES
  time1 = timer_->wallTime();
#endif

  if( !Teuchos::is_null(transform_) )
  {
    if( !tProblem_ )
    {
      tProblem_ = &((*transform_)( *problem_ ));
      tProblem_->SetPDL(unsure);
    }
    swap( tProblem_, problem_ );
    transform_->fwd();
  }


#if defined(Xyce_DETAILED_LINSOLVE_TIMES) || defined(Xyce_VERBOSE_LINEAR)
  time2 = timer_->wallTime();
#endif
#ifdef Xyce_DETAILED_LINSOLVE_TIMES
  cout << "Fwd Trans Time: " << time2-time1 << endl;
#endif

  if( Teuchos::is_null(solver_) ) {

    if (lasProblem_.matrixFree()) {
      Xyce::Report::DevelFatal0().in("N_LAS_ShyLUSolver::solve()") << "cannot work on matrix-free linear systems!";
    }

    Epetra_CrsMatrix * epetraA = dynamic_cast<Epetra_CrsMatrix*>(problem_->GetMatrix());
    solver_ = Teuchos::rcp( new Ifpack_ShyLU( epetraA ) );

    // Set the parameters
    IFPACK_CHK_ERR(solver_->SetParameters(*shyluParams_));

    // Compute symbolic factorization stage of preconditioner.
    IFPACK_CHK_ERR(solver_->Initialize());

#if defined(Xyce_DETAILED_LINSOLVE_TIMES) || defined(Xyce_VERBOSE_LINEAR)
    time1 = timer_->wallTime();
    Xyce::lout() << "ShyLU Initialize Time: " << (time1 - time2) << std::endl;
#endif

    
  }
      
  // Output the linear system to a Matrix Market file every outputLS_ calls if outputLS_ > 0
  static int file_number = 1, base_file_number = 1;
  if (outputLS_ && !lasProblem_.matrixFree()) {
    if (!(file_number % outputLS_)) {
      char file_name[40];
      if (!ReuseFactors) {
        if (file_number == 1) {
          EpetraExt::BlockMapToMatrixMarketFile( "Transformed_BlockMap.mm", (problem_->GetMatrix())->Map() );
        }
        sprintf( file_name, "Transformed_Matrix%d.mm", file_number );
        std::string sandiaReq = "Sandia National Laboratories is a multi-program laboratory managed and operated by Sandia Corporation,\n%";
        sandiaReq += " a wholly owned subsidiary of Lockheed Martin Corporation, for the U.S. Department of Energy's National Nuclear \n%";
        sandiaReq += " Security Administration under contract DE-AC04-94AL85000.\n%\n% Xyce circuit matrix.\n%%";
        EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(problem_->GetMatrix()), sandiaReq.c_str() );
        sprintf( file_name, "Transformed_RHS%d.mm", file_number );
        EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_->GetRHS()) );
      }
    }
    file_number++;
  }
  if (outputBaseLS_ && !lasProblem_.matrixFree()) {
    if (!(base_file_number % outputBaseLS_)) {
      char file_name[40];
      if (!ReuseFactors) {
        if (base_file_number == 1) {
          EpetraExt::BlockMapToMatrixMarketFile( "Base_BlockMap.mm", (tProblem_->GetMatrix())->Map() );
        }
        sprintf( file_name, "Base_Matrix%d.mm", base_file_number );
        std::string sandiaReq = "Sandia National Laboratories is a multi-program laboratory managed and operated by Sandia Corporation,\n%";
        sandiaReq += " a wholly owned subsidiary of Lockheed Martin Corporation, for the U.S. Department of Energy's National Nuclear \n%";
        sandiaReq += " Security Administration under contract DE-AC04-94AL85000.\n%\n% Xyce circuit matrix.\n%%";
        EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(tProblem_->GetMatrix()), sandiaReq.c_str() );
        sprintf( file_name, "Base_RHS%d.mm", base_file_number );
        EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(tProblem_->GetRHS()) );
      }
    }
    base_file_number++;
  }

#if defined(Xyce_DETAILED_LINSOLVE_TIMES) || defined(Xyce_VERBOSE_LINEAR)
  time1 = timer_->wallTime();
#endif

  // Compute the preconditioner
  IFPACK_CHK_ERR(solver_->Compute());
  
#if defined(Xyce_DETAILED_LINSOLVE_TIMES) || defined(Xyce_VERBOSE_LINEAR)
  time2 = timer_->wallTime();
  Xyce::lout() << "ShyLU Compute Time: " << (time2 - time1) << std::endl;
#endif

  // Solve the linear system
  int solveRet = solver_->ApplyInverse( *problem_->GetRHS(), *problem_->GetLHS() );

  if (solveRet)
    Xyce::Report::DevelFatal0().in("N_LAS_ShyLUSolver::solve()") << "ShyLU solver could not be applied!";
  
  //numLinearIters_ = solver_->getNumIters();

#if defined(Xyce_DETAILED_LINSOLVE_TIMES) || defined(Xyce_VERBOSE_LINEAR)
  time1 = timer_->wallTime();
  Xyce::lout() << "ShyLU Solve Time: " << (time1 - time2) << std::endl;
#endif

#ifdef Xyce_DEBUG_LINEAR
  Epetra_MultiVector* b = problem_->GetRHS();
  int numrhs = b->NumVectors();
  std::vector<double> actual_resids( numrhs ), rhs_norm( numrhs );
  Epetra_MultiVector resid( b->Map(), numrhs );
  problem_->GetOperator()->Apply( *(problem_->GetLHS()), resid );
  resid.Update( -1.0, *b, 1.0 );
  resid.Norm2( &actual_resids[0] );
  b->Norm2( &rhs_norm[0] );
  for (int i=0; i<numrhs; i++ ) {
    Xyce::lout() << "Problem " << i << " : \t" <<(actual_resids[i]/rhs_norm[i]) << std::endl;
  }
#endif

  if( !Teuchos::is_null(transform_) )
  {
    transform_->rvs();
    swap( tProblem_, problem_ );
  }

#ifdef Xyce_DETAILED_LINSOLVE_TIMES
  time2 = timer_->wallTime();
  cout << "Rvs Trans Time: " << time2-time1 << endl;
#endif

  // Update the total solution time
  solutionTime_ = timer_->elapsedTime();
#ifdef Xyce_DETAILED_LINSOLVE_TIMES
  cout << "Solve Time: " << solutionTime_ << endl;
#endif

  linearResidual_ = innerTol_/10;
  return 0;
}
