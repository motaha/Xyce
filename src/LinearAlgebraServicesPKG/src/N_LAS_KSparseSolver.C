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
// Filename       : $RCSfile: N_LAS_KSparseSolver.C,v $
//
// Purpose        : KSparse direct solver wrapper
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
// Revision Number: $Revision: 1.17 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <Epetra_LinearProblem.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

// ---------- Xyce Includes ----------

#include <N_LAS_KSparseSolver.h>

#include <N_LAS_Problem.h>

#include <N_LAS_TransformTool.h>

#include <N_UTL_Timer.h>
#include <N_UTL_OptionBlock.h>

#include <N_ERH_ErrorMgr.h>

#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_BlockMapOut.h>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Utils.hpp>

#include <Epetra_CrsKundertSparse.h>
//-----------------------------------------------------------------------------
// Function      : N_LAS_KSparseSolver::N_LAS_KSparseSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
N_LAS_KSparseSolver::N_LAS_KSparseSolver( N_LAS_Problem & prob,
                                          N_UTL_OptionBlock & options )
 : N_LAS_Solver(false),
   lasProblem_(prob),
   problem_(prob.epetraObj()),
   outputLS_(0),
   outputBaseLS_(0),
   outputFailedLS_(0),
   tProblem_(0),
   options_( new N_UTL_OptionBlock( options ) ),
   timer_( new N_UTL_Timer( problem_.GetLHS()->Comm() ) )
{
  setOptions( options );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_KSparseSolver::~N_LAS_KSparseSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
N_LAS_KSparseSolver::~N_LAS_KSparseSolver()
{
  if( timer_ )     delete timer_;
  if( options_ )   delete options_;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_KSparseSolver::setOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_KSparseSolver::setOptions( const N_UTL_OptionBlock & OB )
{
  for( std::list<N_UTL_Param>::const_iterator it_tpL = OB.getParams().begin();
         it_tpL != OB.getParams().end(); ++it_tpL )
  {
    string tag = it_tpL->uTag();

    if( tag == "OUTPUT_LS" ) outputLS_ = it_tpL->getImmutableValue<int>();
    
    if( tag == "OUTPUT_BASE_LS" ) outputBaseLS_ = it_tpL->getImmutableValue<int>();
    
    if( tag == "OUTPUT_FAILED_LS" ) outputFailedLS_ = it_tpL->getImmutableValue<int>();
  }

  if( options_ ) delete options_;
  options_ = new N_UTL_OptionBlock( OB );

#ifdef Xyce_PARALLEL_MPI
  options_->getParams().push_back( N_UTL_Param( "TR_reindex", 1 ) );

  // Turn off partitioning and AMD if we're doing a parallel load serial solve
  options_->getParams().push_back( N_UTL_Param( "TR_partition", 0 ) );
  options_->getParams().push_back( N_UTL_Param( "TR_amd", 0 ) );
#endif

  if( !transform_.get() ) transform_ = N_LAS_TransformTool()( *options_ );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_KSparseSolver::setDefaultOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_KSparseSolver::setDefaultOptions()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_KSparseSolver::setDefaultOption
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_KSparseSolver::setDefaultOption( const string & option )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_KSparseSolver::setParam
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_KSparseSolver::setParam( const N_UTL_Param & param )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_KSparseSolver::getInfo
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_KSparseSolver::getInfo( N_UTL_Param & info )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_KSparseSolver::solve
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
int N_LAS_KSparseSolver::solve( bool ReuseFactors )
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
    file_number++;
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
    base_file_number++;
  }
 
#ifdef Xyce_DEBUG_LINEAR
    // Set the traceback mode in Epetra so it prints out warnings
    dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix())->SetTracebackMode( 2 );
#else
    dynamic_cast<Epetra_CrsMatrix*>(prob->GetMatrix())->SetTracebackMode( 0 );
#endif

    // Import the parallel matrix to a serial one, if necessary.
    prob = importToSerial();

    // Create solver if one doesn't exist.
    if (solver_ == Teuchos::null)
      solver_ = Teuchos::rcp( new Epetra_CrsKundertSparse( prob ) );

  // Perform linear solve using factorization
#ifdef Xyce_VERBOSE_LINEAR
  double begSolveTime = timer_->elapsedTime();
#endif      

    linearStatus = solver_->Solve( !ReuseFactors );
 
    // Export solution back to global system, if necessary.
    exportToGlobal();

#ifdef Xyce_VERBOSE_LINEAR
    double endSolveTime = timer_->elapsedTime();
    Xyce::lout() << "  KSparse Solve Time: " << (endSolveTime - begSolveTime) << std::endl;
#endif      

    if (linearStatus != 0) {

      // Inform user that singular matrix was found and linear solve has failed.
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_WARNING_0,
        "Numerically singular matrix found by KSparse, returning zero solution to nonlinear solver!");

      // Put zeros in the solution since KSparse was not able to solve this problem
      prob->GetLHS()->PutScalar( 0.0 );
      // Output the singular linear system to a Matrix Market file if outputFailedLS_ > 0
      if (outputFailedLS_) {
        failure_number++;
        char file_name[40];
        if (failure_number == 1) {
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
    }

#ifdef Xyce_DEBUG_LINEAR
    double resNorm = 0.0, bNorm = 0.0;
    Epetra_MultiVector res( prob->GetLHS()->Map(), prob->GetLHS()->NumVectors() );
    prob->GetOperator()->Apply( *(prob->GetLHS()), res );
    res.Update( 1.0, *(prob->GetRHS()), -1.0 );
    res.Norm2( &resNorm );
    prob->GetRHS()->Norm2( &bNorm );
    Xyce::lout() << "Linear System Residual (KSparse) : " << (resNorm/bNorm) << std::endl;
#endif

 if( transform_.get() ) transform_->rvs();
  
  // Update the total solution time
  solutionTime_ = timer_->elapsedTime();

#ifdef Xyce_VERBOSE_LINEAR
  Xyce::lout() << "Total Linear Solution Time (KSparse): " << solutionTime_ << std::endl;
#endif

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_KSparseSolver::importToSerial
// Purpose       : Import a distributed matrix to a serial one for the solver
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystems Modeling
// Creation Date : 07/01/10
//-----------------------------------------------------------------------------
Epetra_LinearProblem * N_LAS_KSparseSolver::importToSerial()
{
#ifdef Xyce_PARALLEL_MPI
    Epetra_CrsMatrix * origMat = dynamic_cast<Epetra_CrsMatrix *>(problem_.GetOperator());

    // If we haven't set up the maps for serial matrix storage on proc 0, then do it now.
    if (serialMap_ == Teuchos::null) {
      const Epetra_Map& origMap = origMat->RowMap();
      int MyPID = origMap.Comm().MyPID(); 
      int NumGlobalElements = origMap.NumGlobalElements();
      int NumMyElements = NumGlobalElements;
      if (MyPID != 0)
        NumMyElements = 0;
      serialMap_ = Teuchos::rcp(new Epetra_Map(-1, NumMyElements, 0, origMap.Comm()));
      int NumVectors = problem_.GetRHS()->NumVectors() ;
    
      serialLHS_ = Teuchos::rcp( new Epetra_MultiVector( *serialMap_, NumVectors ));
      serialRHS_ = Teuchos::rcp (new Epetra_MultiVector( *serialMap_, NumVectors ));
      serialImporter_ = Teuchos::rcp(new Epetra_Import( *serialMap_, origMap ));
      serialMat_ = Teuchos::rcp( new Epetra_CrsMatrix( Copy, *serialMap_, 0 )) ;
      serialProblem_ = Teuchos::rcp( new Epetra_LinearProblem( &*serialMat_, &*serialLHS_, &*serialRHS_ ) );
    }

    // Import linear system from problem
    serialMat_->Import( *origMat, *serialImporter_, Insert );
    serialMat_->FillComplete();
    serialMat_->OptimizeStorage(); 
    serialLHS_->Import( *problem_.GetLHS(), *serialImporter_, Insert );
    serialRHS_->Import( *problem_.GetRHS(), *serialImporter_, Insert );   

    return &*serialProblem_;
#else
    return &problem_;
#endif
} 

//-----------------------------------------------------------------------------
// Function      : N_LAS_KSparseSolver::exportToGlobal
// Purpose       : Export the serial solution to the global one
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystems Modeling
// Creation Date : 07/01/10
//-----------------------------------------------------------------------------
int N_LAS_KSparseSolver::exportToGlobal()
{
#ifdef Xyce_PARALLEL_MPI
  // Return solution back to global problem
  problem_.GetLHS()->Export( *serialLHS_, *serialImporter_, Insert ) ;
#endif
  return 0;
}
