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
// Filename       : $RCSfile: N_LAS_BelosSolver.C,v $
//
// Purpose        : Implementation file for the Belos linear solver interface.
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
// Revision Number: $Revision: 1.39 $
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

#include <N_UTL_fwd.h>
#include <N_UTL_Misc.h>

#include <N_LAS_BelosSolver.h>

#include <N_UTL_OptionBlock.h>

#include <N_UTL_Timer.h>

#include <N_ERH_ErrorMgr.h>

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_InvOperator.h>

#include <BelosBlockGmresSolMgr.hpp>
#include <BelosGCRODRSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosLinearProblem.hpp>

#include <N_LAS_TransformTool.h>
#include <N_LAS_Problem.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Preconditioner.h>
#include <N_LAS_TrilinosPrecondFactory.h>

#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_BlockMapOut.h>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Utils.hpp>

//-----------------------------------------------------------------------------
// Function      : N_LAS_BelosSolver::~N_LAS_BelosSolver
// Purpose       : Default destructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
N_LAS_BelosSolver::~N_LAS_BelosSolver()
{
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BelosSolver::N_LAS_BelosSolver
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
N_LAS_BelosSolver::N_LAS_BelosSolver( N_LAS_Problem & problem,
				      N_UTL_OptionBlock & options )
: N_LAS_Solver(true),
#ifdef Xyce_VERBOSE_LINEAR
  outputDefault_(Belos::Errors + Belos::Warnings + Belos::StatusTestDetails),
#else
  outputDefault_(0),
#endif
  maxIterDefault_(500),
  KSpaceDefault_(500),
  toleranceDefault_(1.0e-12),
  recycleDefault_(10),
  belosSolverDefault_("Block GMRES"),
  outputLS_(0),
  outputBaseLS_(0),
  lasProblem_(problem),
  problem_(&(problem.epetraObj())),
  updatedParams_(false),
  linearResidual_(1.0),
  tProblem_(0),
  isPrecSet_(false)
{
  options_ = Teuchos::rcp( new N_UTL_OptionBlock( options ) );
  timer_ = Teuchos::rcp( new N_UTL_Timer( problem_->GetLHS()->Comm() ) );

  belosParams_ = Teuchos::rcp( new Teuchos::ParameterList() );
  setDefaultOptions();
  setOptions( *options_ );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BelosSolver::setDefaultOptions
// Purpose       : resets Aztec options
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 9/25/07
//-----------------------------------------------------------------------------
bool N_LAS_BelosSolver::setDefaultOptions()
{
  output_ = outputDefault_;
  maxIter_ = maxIterDefault_;
  numLinearIters_ = 0;
  KSpace_ = KSpaceDefault_;
  recycle_ = recycleDefault_;
  tolerance_ = toleranceDefault_;
  belosSolver_ = belosSolverDefault_;
  updatedParams_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BelosSolver::setDefaultOption
// Purpose       : resets Aztec option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 9/25/07
//-----------------------------------------------------------------------------
bool N_LAS_BelosSolver::setDefaultOption( const std::string & option )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BelosSolver::setOptions
// Purpose       : sets Aztec options and params from modelblock
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool N_LAS_BelosSolver::setOptions(const N_UTL_OptionBlock& OB)
{
  belosParams_->set("Verbosity",output_);
#ifdef Xyce_VERBOSE_LINEAR
  belosParams_->set("Output Frequency", 50);
  belosParams_->set("Output Style", Belos::Brief);
#endif
  belosParams_->set("Maximum Iterations",maxIter_);
  belosParams_->set("Num Blocks", KSpace_);
  belosParams_->set("Block Size", 1);
  belosParams_->set("Convergence Tolerance", tolerance_);
  belosParams_->set("Orthogonalization", "ICGS");

  // Number of vectors in recycle space for GCRODR
  belosParams_->set("Num Recycled Blocks", recycle_);

  std::list<N_UTL_Param>::const_iterator it_tpL = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator end_tpL = OB.getParams().end();
  for (; it_tpL != end_tpL; ++it_tpL)
  {
    setParam( *it_tpL );
  }

  // store for restart of solver_
  if( &OB != &*options_ )
  {
    options_ = Teuchos::rcp( new N_UTL_OptionBlock(OB) );
  }

  // set singleton filtering as default for iterative solvers
  if (!OB.tagExists("TR_singleton_filter"))
  {
    options_->getParams().push_back( N_UTL_Param( "TR_singleton_filter", 1 ) );
  }

  if (!lasProblem_.matrixFree())
  {
    // create the transformation object if needed.
    if( Teuchos::is_null(transform_) ) transform_ = N_LAS_TransformTool()( OB );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BelosSolver::setParam
// Purpose       : sets Aztec option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool N_LAS_BelosSolver::setParam( const N_UTL_Param & param )
{
  std::string tag = param.tag();
  std::string uTag = param.uTag();

  // Set our copies of these parameters.
  if( tag == "AZ_max_iter" )
    setMaxIter(param.getImmutableValue<int>());
  if( tag == "AZ_kspace" )
    setKSpace(param.getImmutableValue<int>());
  else if( tag == "AZ_tol" )
    setTolerance(param.getImmutableValue<double>());
  else if( uTag == "OUTPUT_LS" )
    outputLS_ = param.getImmutableValue<int>();
  else if( uTag == "OUTPUT_BASE_LS" )
    outputBaseLS_ = param.getImmutableValue<int>();
  else if( tag == "BELOS_SOLVER_TYPE" )
   belosSolver_ = param.usVal();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BelosSolver::getInfo
// Purpose       : sets Aztec option
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool N_LAS_BelosSolver::getInfo( N_UTL_Param & param )
{
  if( param.tag() == "AZ_max_iter" )
    param.setVal( maxIter_ );
  else if( param.tag() == "Iterations" )
    param.setVal( (int)numLinearIters_ );
  else if( param.tag() == "AZ_tol" )
    param.setVal( tolerance_ );
  else
    return false;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BelosSolver::setBelosOption
// Purpose       : sets Belos option
// Special Notes : Takes a string as the option identifier
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Computational Sciences
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool N_LAS_BelosSolver::setBelosOption_(const char * paramName,
                                            const int val)
{
  return setBelosCntl_( N_UTL_Param(paramName, val) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BelosSolver::setBelosParam
// Purpose       : sets Belos parameter
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool N_LAS_BelosSolver::setBelosParam_(const char * paramName,
                                         const double val)
{
  return setBelosCntl_( N_UTL_Param(paramName, val) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Belos::printParams_
// Purpose       : Print out the linear solver parameter values.
// Special Notes :
// Scope         : Private
// Creator       : Heidi Thornquist, SNL, Computational Sciences
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
void N_LAS_BelosSolver::printParams_() const
{}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BelosSolver::solve
// Purpose       : Calls the actual solver to solve Ax=b.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
int N_LAS_BelosSolver::solve( bool ReuseFactors )
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
    std::swap( tProblem_, problem_ );
    transform_->fwd();
  }


#ifdef Xyce_DETAILED_LINSOLVE_TIMES
  time2 = timer_->wallTime();
  Xyce::dout() << "Fwd Trans Time: " << time2-time1 << std::endl;
#endif

  if( Teuchos::is_null(solver_) ) {

    Teuchos::RCP< Epetra_Operator > A;
    if (lasProblem_.matrixFree()) {
      A = Teuchos::rcp( problem_->GetOperator(), false );
    }
    else {
      A = Teuchos::rcp( problem_->GetMatrix(), false );
    }
    Teuchos::RCP< Epetra_MultiVector> X = Teuchos::rcp( problem_->GetLHS(), false );
    Teuchos::RCP< Epetra_MultiVector> B = Teuchos::rcp( problem_->GetRHS(), false );
    belosProblem_ =
      Teuchos::rcp( new Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>( A, X, B ) );

    // Create Belos solver (Block GMRES is default)
    if (belosSolver_ == "GCRODR") {
      solver_ = Teuchos::rcp( new Belos::GCRODRSolMgr<double,Epetra_MultiVector,Epetra_Operator>() );
    }
    else {
      solver_ = Teuchos::rcp( new Belos::BlockGmresSolMgr<double,Epetra_MultiVector,Epetra_Operator>() );
    }

    // Reduce the size of the Krylov space if the linear system is smaller than the default size
    // to avoid warning messages.
    int lsDim = X->GlobalLength();
    if (lsDim < KSpace_) {
      KSpace_ = lsDim;
      belosParams_->set("Num Blocks", KSpace_);
    }
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
    // file_number++;  This is incremented below after the solution is written to file.
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
    // base_file_number++;  This is incremented below after the solution is written to file.
  }

#ifdef Xyce_DETAILED_LINSOLVE_TIMES
  time1 = timer_->wallTime();
#endif

  Teuchos::RCP<N_LAS_Problem> tmpProblem = Teuchos::rcp( new N_LAS_Problem( Teuchos::rcp(problem_,false) ) );

  // Create the preconditioner if we don't have one.
  if ( Teuchos::is_null( precond_ ) ) {
    N_LAS_TrilinosPrecondFactory factory( *options_ );
    precond_ = factory.create( tmpProblem );
    isPrecSet_ = false;
  }

  // Initialize the values compute the preconditioner.
  bool initRet = precond_->initValues( tmpProblem );
  if (!initRet)
    Xyce::Report::DevelFatal0().in("N_LAS_Preconditioner::initValues()") << "preconditioner could not be initialized!";

  // Compute the preconditioner
  bool compRet = precond_->compute();
  if (!compRet)
    Xyce::Report::DevelFatal0().in("N_LAS_Preconditioner::compute()") << "preconditioner could not be computed!";

  // Set the preconditioner as an operator if it was just constructed, else inform the Belos
  // solver to use no preconditioning.
  if ( !isPrecSet_ ) {
    if ( !Teuchos::is_null( precond_->epetraObj() ) ) {
      belosPrecond_ = Teuchos::rcp( new Epetra_InvOperator( &*precond_->epetraObj() ) );
    }
      else {
      belosPrecond_ = Teuchos::null;
    }
    belosProblem_->setRightPrec( belosPrecond_ );
    belosProblem_->setProblem();
    solver_->setProblem( belosProblem_ );
    solver_->setParameters( belosParams_ );
    isPrecSet_ = true;
  }

#ifdef Xyce_DETAILED_LINSOLVE_TIMES
  time2 = timer_->wallTime();
  Xyce::dout() << "Ifpack Time: " << time2-time1 << std::endl;
#elif defined(Xyce_VERBOSE_LINEAR)
  time2 = timer_->wallTime();
#endif


  // Reset the problem for the next solve.
  belosProblem_->setProblem();

  // Solve the problem.
  Belos::ReturnType linearStatus;
  linearStatus = solver_->solve();
  numLinearIters_ = solver_->getNumIters();

#if defined(Xyce_DETAILED_LINSOLVE_TIMES) || defined(Xyce_VERBOSE_LINEAR)
  time1 = timer_->wallTime();
  Xyce::lout() << "Belos Solve Time: " << (time1 - time2) << std::endl;
#endif

/*
  Epetra_MultiVector* b = problem_->GetRHS();
  int numrhs = b->NumVectors();
  std::vector<double> actual_resids( numrhs ), rhs_norm( numrhs );
  Epetra_MultiVector resid( b->Map(), numrhs );
  problem_->GetOperator()->Apply( *(problem_->GetLHS()), resid );
  resid.Update( -1.0, *b, 1.0 );
  resid.Norm2( &actual_resids[0] );
  b->Norm2( &rhs_norm[0] );
  for (int i=0; i<numrhs; i++ ) {
    Xyce::dout() << "Problem " << i << " : \t" << actual_resids[i]/rhs_norm[i] << std::endl;
  }
*/

  if( !Teuchos::is_null(transform_) )
  {
    transform_->rvs();
    std::swap( tProblem_, problem_ );
  }

/*
  Epetra_MultiVector* b2 = problem_->GetRHS();
  Epetra_MultiVector resid2( b2->Map(), numrhs );
  problem_->GetOperator()->Apply( *(problem_->GetLHS()), resid2 );
  resid2.Update( -1.0, *b2, 1.0 );
  resid2.Norm2( &actual_resids[0] );
  b2->Norm2( &rhs_norm[0] );
  for (int i=0; i<numrhs; i++ ) {
    Xyce::dout() << "Problem " << i << " : \t" << actual_resids[i]/rhs_norm[i] << std::endl;
  }
*/

#ifdef Xyce_DETAILED_LINSOLVE_TIMES
  time2 = timer_->wallTime();
  Xyce::dout() << "Rvs Trans Time: " << time2-time1 << std::endl;
#endif

  // Output computed solution vectors, if requested.
  if (outputLS_ && !lasProblem_.matrixFree()) {
    if (!(file_number % outputLS_)) {
      char file_name[40];
      sprintf( file_name, "Transformed_Soln%d.mm", file_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(tProblem_->GetLHS()) );
    }
    file_number++;
  }
  if (outputBaseLS_ && !lasProblem_.matrixFree()) {
    if (!(base_file_number % outputBaseLS_)) {
      char file_name[40];
      sprintf( file_name, "Base_Soln%d.mm", base_file_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_->GetLHS()) );
    }
    base_file_number++;
  }

  // Update the total solution time
  solutionTime_ = timer_->elapsedTime();
#ifdef Xyce_DETAILED_LINSOLVE_TIMES
  Xyce::dout() << "Solve Time: " << solutionTime_ << std::endl;
#endif

  // Belos does not return the residual, so use a fake residual that indicates convergence
  // if the solver converged.
  if (linearStatus == Belos::Unconverged)
  {
#ifdef Xyce_VERBOSE_LINEAR
    if (solver_->isLOADetected())
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0,
			   "Belos::BlockGmresSolMgr has detected a loss of accuracy!\n");
    }
#endif
  }

  linearResidual_ = tolerance_/10;
  return 0;
}
