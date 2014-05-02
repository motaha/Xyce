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
// Filename       : $RCSfile: N_LAS_AmesosSolver.C,v $
//
// Purpose        : Amesos direct solver wrapper
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/20/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.58 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <Amesos.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Export.h>
#include <Epetra_Map.h>

// ---------- Xyce Includes ----------

#include <N_LAS_AmesosSolver.h>

#include <N_LAS_Problem.h>

#include <N_LAS_TransformTool.h>

#include <N_UTL_fwd.h>
#include <N_UTL_Timer.h>
#include <N_UTL_OptionBlock.h>

#include <N_ERH_ErrorMgr.h>

#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_BlockMapOut.h>

#include <Teuchos_Utils.hpp>
//-----------------------------------------------------------------------------
// Function      : N_LAS_AmesosSolver::N_LAS_AmesosSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
N_LAS_AmesosSolver::N_LAS_AmesosSolver( const std::string& type,
                                        N_LAS_Problem & prob,
                                        N_UTL_OptionBlock & options )
 : N_LAS_Solver(false),
   type_(type),
   lasProblem_(prob),
   problem_(prob.epetraObj()),
   solver_(0),
   repivot_(true),
   reindex_(false),
   outputLS_(0),
   outputBaseLS_(0),
   outputFailedLS_(0),
   tProblem_(0),
   optProb_(0),
   optMat_(0),
   origMat_(0),
   optExporter_(0),
   options_( new N_UTL_OptionBlock( options ) ),
   timer_( new N_UTL_Timer( problem_.GetLHS()->Comm() ) )
{
  setOptions( options );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_AmesosSolver::~N_LAS_AmesosSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
N_LAS_AmesosSolver::~N_LAS_AmesosSolver()
{
  if(solver_)    delete solver_;
  if(timer_)     delete timer_;
  if(options_)   delete options_;
  if( optProb_ ) delete optProb_;
  if( optMat_ )  delete optMat_;
  if( optExporter_ ) delete optExporter_;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_AmesosSolver::setOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_AmesosSolver::setOptions( const N_UTL_OptionBlock & OB )
{
  for( std::list<N_UTL_Param>::const_iterator it_tpL = OB.getParams().begin();
         it_tpL != OB.getParams().end(); ++it_tpL )
  {
    std::string tag = it_tpL->uTag();

    if( tag == "KLU_REPIVOT" ) repivot_ = static_cast<bool>(it_tpL->getImmutableValue<int>());
    
    if( tag == "KLU_REINDEX" ) reindex_ = static_cast<bool>(it_tpL->getImmutableValue<int>());

    if( tag == "OUTPUT_LS" ) outputLS_ = it_tpL->getImmutableValue<int>();

    if( tag == "OUTPUT_BASE_LS" ) outputBaseLS_ = it_tpL->getImmutableValue<int>();

    if( tag == "OUTPUT_FAILED_LS" ) outputFailedLS_ = it_tpL->getImmutableValue<int>();
  }

  if( options_ ) delete options_;
  options_ = new N_UTL_OptionBlock( OB );

#ifdef Xyce_PARALLEL_MPI
  options_->getParams().push_back( N_UTL_Param( "TR_reindex", 1 ) );

  // Turn off partitioning and AMD if we're doing a parallel load serial solve
  if (type_ == "KLU") {
    options_->getParams().push_back( N_UTL_Param( "TR_partition", 0 ) );
    options_->getParams().push_back( N_UTL_Param( "TR_amd", 0 ) );
  }
#endif

  if( !transform_.get() ) transform_ = N_LAS_TransformTool()( *options_ );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_AmesosSolver::setDefaultOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_AmesosSolver::setDefaultOptions()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_AmesosSolver::setDefaultOption
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_AmesosSolver::setDefaultOption( const std::string & option )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_AmesosSolver::setParam
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_AmesosSolver::setParam( const N_UTL_Param & param )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_AmesosSolver::getInfo
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_AmesosSolver::getInfo( N_UTL_Param & info )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_AmesosSolver::solve
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
int N_LAS_AmesosSolver::solve( bool ReuseFactors )
{
  // Start the timer...
  timer_->resetStartTime();

  int linearStatus = 0;

  Epetra_LinearProblem * prob = &problem_;

  if( transform_.get() )
  {
    if( !tProblem_ )
      tProblem_ = &((*transform_)( problem_ ));
    prob = tProblem_;
    transform_->fwd();
  }

  // Output the linear system to a Matrix Market file every outputLS_ calls if outputLS_ > 0
  static int failure_number = 0, file_number = 1, base_file_number = 1;
  if (outputLS_) {
    if (!(file_number % outputLS_)) {
      char file_name[40];
      if (!ReuseFactors) {
        if (file_number == 1) {
          EpetraExt::BlockMapToMatrixMarketFile( "Transformed_BlockMap.mm", (prob->GetMatrix())->Map() );
        }
        sprintf( file_name, "Transformed_Matrix%d.mm", file_number );

        std::string sandiaReq = "Sandia National Laboratories is a multi-program laboratory managed and operated by Sandia Corporation,\n%";
        sandiaReq += " a wholly owned subsidiary of Lockheed Martin Corporation, for the U.S. Department of Energy's National Nuclear \n%";
        sandiaReq += " Security Administration under contract DE-AC04-94AL85000.\n%\n% Xyce circuit matrix.\n%%";

        EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(prob->GetMatrix()), sandiaReq.c_str() );
        sprintf( file_name, "Transformed_RHS%d.mm", file_number );
        EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(prob->GetRHS()) );
      }
    }
    // file_number++;  This will be incremented after the solution vector is written to file.
  }
  if (outputBaseLS_) {
    if (!(base_file_number % outputBaseLS_)) {
      char file_name[40];
      if (!ReuseFactors) {
        if (base_file_number == 1) {
          EpetraExt::BlockMapToMatrixMarketFile( "Base_BlockMap.mm", (problem_.GetMatrix())->Map() );
        }
        sprintf( file_name, "Base_Matrix%d.mm", base_file_number );

        std::string sandiaReq = "Sandia National Laboratories is a multi-program laboratory managed and operated by Sandia Corporation,\n%";
        sandiaReq += " a wholly owned subsidiary of Lockheed Martin Corporation, for the U.S. Department of Energy's National Nuclear \n%";
        sandiaReq += " Security Administration under contract DE-AC04-94AL85000.\n%\n% Xyce circuit matrix.\n%%";

        EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(problem_.GetMatrix()), sandiaReq.c_str() );
        sprintf( file_name, "Base_RHS%d.mm", base_file_number );
        EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_.GetRHS()) );
      }
    }
    // base_file_number++;  This will be incremented after the solution vector is written to file.
  }

#ifdef Xyce_DEBUG_LINEAR
    // Set the traceback mode in Epetra so it prints out warnings
    dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix())->SetTracebackMode( 2 );
#else
    dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix())->SetTracebackMode( 0 );
#endif

  Amesos localAmesosObject;
  if( !solver_ )
  {
    // the Query() function expects a string
    // in lower case with the first letter in upper case
    // So, our "KLU" must become "Klu"
    std::string solverType( type_ );
    if( type_ == "KLU" )
    {
        solverType = "Amesos_Klu";
    }
    else if( type_ == "SUPERLU" )
    {
	solverType = "Amesos_Superlu";
    }
    else if( type_ == "SUPERLUDIST" )
    {
	solverType = "Amesos_Superludist";
    }
    else if( type_ == "PARAKLETE" )
    {
        solverType = "Amesos_Paraklete";
    }
    else if( type_ == "PARDISO" )
    {
        solverType = "Amesos_Pardiso";
    }
    else if( type_ == "LAPACK" )
    {
        solverType = "Amesos_Lapack";
    }

    if( !localAmesosObject.Query( solverType ) )
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,
        "Unknown or Unavailable Linear Solver: " + type_ );


#ifndef Xyce_PARALLEL_MPI
    //setup optimized storage version of problem for serial
    //only do this if the linear system is nontrivial (not a single equation)
    origMat_ = dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix());
    if (origMat_->NumGlobalRows() > 1) {
      Epetra_Map const& rowMap = origMat_->RowMap();
      Epetra_BlockMap const& blockRowMap = dynamic_cast<Epetra_BlockMap const&>(rowMap);
      optMat_ = new Epetra_CrsMatrix( Copy, rowMap, 0 );
      optExporter_ = new Epetra_Export( blockRowMap, blockRowMap );
      optMat_->Export( *origMat_, *optExporter_, Insert );
      optMat_->FillComplete();
      optMat_->OptimizeStorage();

      optProb_ = new Epetra_LinearProblem( optMat_, prob->GetLHS(), prob->GetRHS() );
      prob = optProb_;
    }
#endif

    solver_ = localAmesosObject.Create( solverType, *prob );

    Teuchos::ParameterList params;

#ifndef Xyce_PARALLEL_MPI
    // Inform solver not to check inputs to reduce overhead.
    params.set( "TrustMe", true );
    // If repivot == true (default), recompute the pivot order each numeric factorization,
    // else try to re-use pivot order to expedite numeric factorization.
    params.set( "Refactorize", !repivot_ );
#else
    if (type_ == "SUPERLUDIST") {
      Teuchos::ParameterList& sludistParams = params.sublist("Superludist");
      sludistParams.set("ReuseSymbolic", true );
    }
#endif

    // Let Amesos reindex the linear problem.
    // NOTE:  This is used by MPDE and HB since the map indices are not continguous.
    if (reindex_) {
      params.set( "Reindex", reindex_ );
    }

#ifdef Xyce_VERBOSE_LINEAR
    Xyce::lout() << "N_LAS_AmesosSolver::solve() setting solver : " << type_ << "\n"
                 << "N_LAS_AmesosSolver::solve() setting parameters : " << params << std::endl;
#endif

    solver_->SetParameters( params );

#ifdef Xyce_VERBOSE_LINEAR
    double begSymTime = timer_->elapsedTime();
#endif

    // Perform symbolic factorization and check return value for failure
    linearStatus = solver_->SymbolicFactorization();
    if (linearStatus != 0) return linearStatus;

#ifdef Xyce_VERBOSE_LINEAR
    double endSymTime = timer_->elapsedTime();
    Xyce::lout() << "  Amesos (" << type_ << ") Symbolic Factorization Time: "
                 << (endSymTime - begSymTime) << std::endl;
#endif
  }

  if( optMat_ ) optMat_->Export( *origMat_, *optExporter_, Insert );

  // Perform numeric factorization and check return value for failure
  if( !ReuseFactors ) {
#ifdef Xyce_VERBOSE_LINEAR
    double begNumTime = timer_->elapsedTime();
#endif
    linearStatus = solver_->NumericFactorization();
#ifdef Xyce_VERBOSE_LINEAR
    double endNumTime = timer_->elapsedTime();
    Xyce::lout() << "  Amesos (" << type_ << ") Numeric Factorization Time: "
                 << (endNumTime - begNumTime) << std::endl;
#endif
    if (linearStatus != 0) {

      // Inform user that singular matrix was found and linear solve has failed.
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_WARNING_0,
        "Numerically singular matrix found by Amesos, returning zero solution to nonlinear solver!");

      // Put zeros in the solution since Amesos was not able to solve this problem
      prob->GetLHS()->PutScalar( 0.0 );
      // Output the singular linear system to a Matrix Market file if outputFailedLS_ > 0
      if (outputFailedLS_) {
        failure_number++;
        char file_name[40];
        if (failure_number== 1) {
          EpetraExt::BlockMapToMatrixMarketFile( "Failed_BlockMap.mm", (prob->GetMatrix())->Map() );
        }
        sprintf( file_name, "Failed_Matrix%d.mm", failure_number );
        std::string sandiaReq = "Sandia National Laboratories is a multi-program laboratory managed and operated by Sandia Corporation,\n%";
        sandiaReq += " a wholly owned subsidiary of Lockheed Martin Corporation, for the U.S. Department of Energy's National Nuclear \n%";
        sandiaReq += " Security Administration under contract DE-AC04-94AL85000.\n%\n% Xyce circuit matrix.\n%%";

        EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(prob->GetMatrix()), sandiaReq.c_str() );
        sprintf( file_name, "Failed_RHS%d.mm", failure_number );
        EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(prob->GetRHS()) );
      }
      return 0;  // return 0 instead of linearStatus and let nonlinear solver decide what to do.
    }
  }

  // Perform linear solve using factorization
#ifdef Xyce_VERBOSE_LINEAR
  double begSolveTime = timer_->elapsedTime();
#endif

  solver_->Solve();

#ifdef Xyce_VERBOSE_LINEAR
    double endSolveTime = timer_->elapsedTime();
    Xyce::lout() << "  Amesos (" << type_ << ") Solve Time: "
                 << (endSolveTime - begSolveTime) << std::endl;
#endif

#ifdef Xyce_DEBUG_LINEAR
    double resNorm = 0.0, bNorm = 0.0;
    Epetra_MultiVector res( prob->GetLHS()->Map(), prob->GetLHS()->NumVectors() );
    prob->GetOperator()->Apply( *(prob->GetLHS()), res );
    res.Update( 1.0, *(prob->GetRHS()), -1.0 );
    res.Norm2( &resNorm );
    prob->GetRHS()->Norm2( &bNorm );
    Xyce::lout() << "Linear System Residual (AMESOS_" << type_ << "): "
                 << (resNorm/bNorm) << std::endl;
#endif

  if( transform_.get() ) transform_->rvs();

  // Output computed solution vectors, if requested.
  if (outputLS_) {
    if (!(file_number % outputLS_)) {
      char file_name[40];
      sprintf( file_name, "Transformed_Soln%d.mm", file_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_.GetLHS()) );
    }
    file_number++;
  }
  if (outputBaseLS_) {
    if (!(base_file_number % outputBaseLS_)) {
      char file_name[40];
      sprintf( file_name, "Base_Soln%d.mm", base_file_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(prob->GetLHS()) );
    }
    base_file_number++;
  }

  // Update the total solution time
  solutionTime_ = timer_->elapsedTime();

#ifdef Xyce_VERBOSE_LINEAR
  Xyce::lout() << "Total Linear Solution Time (Amesos " << type_ << "): "
               << solutionTime_ << std::endl;
#endif

  return 0;
}
