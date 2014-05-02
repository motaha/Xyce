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
// Filename       : $RCSfile: N_LAS_SimpleSolver.C,v $
//
// Purpose        : Simple direct solver when the matrix is 1x1
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 03/07/2013
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

#include <Epetra_LinearProblem.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Export.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

// ---------- Xyce Includes ----------

#include <N_LAS_SimpleSolver.h>

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
// Function      : N_LAS_SimpleSolver::N_LAS_SimpleSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
N_LAS_SimpleSolver::N_LAS_SimpleSolver( N_LAS_Problem & prob,
                                        N_UTL_OptionBlock & options )
 : N_LAS_Solver(false),
   lasProblem_(prob),
   problem_(prob.epetraObj()),
   outputLS_(0),
   outputBaseLS_(0),
   outputFailedLS_(0),
   options_( new N_UTL_OptionBlock( options ) ),
   timer_( new N_UTL_Timer( problem_.GetLHS()->Comm() ) )
{
  setOptions( options );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_SimpleSolver::~N_LAS_SimpleSolver
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
N_LAS_SimpleSolver::~N_LAS_SimpleSolver()
{
  if(timer_)     delete timer_;
  if(options_)   delete options_;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_SimpleSolver::setOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_SimpleSolver::setOptions( const N_UTL_OptionBlock & OB )
{
  for( std::list<N_UTL_Param>::const_iterator it_tpL = OB.getParams().begin();
         it_tpL != OB.getParams().end(); ++it_tpL )
  {
    std::string tag = it_tpL->uTag();

    if( tag == "OUTPUT_LS" ) outputLS_ = it_tpL->getImmutableValue<int>();

    if( tag == "OUTPUT_BASE_LS" ) outputBaseLS_ = it_tpL->getImmutableValue<int>();

    if( tag == "OUTPUT_FAILED_LS" ) outputFailedLS_ = it_tpL->getImmutableValue<int>();
  }

  if( options_ ) delete options_;
  options_ = new N_UTL_OptionBlock( OB );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_SimpleSolver::setDefaultOptions
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_SimpleSolver::setDefaultOptions()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_SimpleSolver::setDefaultOption
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_SimpleSolver::setDefaultOption( const std::string & option )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_SimpleSolver::setParam
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_SimpleSolver::setParam( const N_UTL_Param & param )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_SimpleSolver::getInfo
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
bool N_LAS_SimpleSolver::getInfo( N_UTL_Param & info )
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_SimpleSolver::solve
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
int N_LAS_SimpleSolver::solve( bool ReuseFactors )
{
  // Start the timer...
  timer_->resetStartTime();

  Epetra_LinearProblem * prob = &problem_;

  // Output the linear system to a Matrix Market file every outputLS_ calls if outputLS_ > 0 (or outputBaseLS_)
  static int failure_number = 0, file_number = 1, base_file_number = 1;
  if (outputBaseLS_ || outputLS_) {
    if (!(base_file_number % outputBaseLS_) || !(file_number % outputLS_)) {
      char file_name[40];
      if (!ReuseFactors) {
        if (base_file_number == 1) {
          EpetraExt::BlockMapToMatrixMarketFile( "Simple_BlockMap.mm", (problem_.GetMatrix())->Map() );
        }
        sprintf( file_name, "Simple_Matrix%d.mm", base_file_number );

        std::string sandiaReq = "Sandia National Laboratories is a multi-program laboratory managed and operated by Sandia Corporation,\n%";
        sandiaReq += " a wholly owned subsidiary of Lockheed Martin Corporation, for the U.S. Department of Energy's National Nuclear \n%";
        sandiaReq += " Security Administration under contract DE-AC04-94AL85000.\n%\n% Xyce circuit matrix.\n%%";

        EpetraExt::RowMatrixToMatrixMarketFile( file_name, *(problem_.GetMatrix()), sandiaReq.c_str() );
        sprintf( file_name, "Simple_RHS%d.mm", base_file_number );
        EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_.GetRHS()) );
      }
    }
    // base_file_number++;  This will be incremented after the solution vector is written to file.
  }

  // Perform linear solve using factorization
#ifdef Xyce_VERBOSE_LINEAR
  double begSolveTime = timer_->elapsedTime();
#endif

  // Make sure this is a trivial linear system.
  int NumGlobalRows = problem_.GetMatrix()->NumGlobalRows();
  if (NumGlobalRows > 1) {
    // Inform user that a nontrivial matrix was found and linear solve has failed.
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_ERROR_0,
                            "Nontrivial matrix has been found, this cannot be handled by this linear solver!");
  }

  int NumIndices = 0;
  int Length = problem_.GetMatrix()->MaxNumEntries();
  std::vector<int> Indices(Length);
  std::vector<double> Values(Length);
  int NumMyRows = problem_.GetMatrix()->NumMyRows();
  int localSingularMat = 0, singularMat = 0;

  for (int i=0; i<NumMyRows; ++i) {

    // Get ith row
    EPETRA_CHK_ERR(problem_.GetMatrix()->ExtractMyRowCopy(i, Length, NumIndices, &Values[0], &Indices[0]));
    if (NumIndices != 1) {
      // Inform user that an empty matrix was found and linear solve has failed.
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_ERROR_0,
                              "Empty matrix has been found, this linear solve has failed!");
    }

    double pivot = Values[0];
    if (pivot != 0.0)
      problem_.GetLHS()->Scale( 1.0/pivot, *(problem_.GetRHS()) );
    else
      localSingularMat = true;
  }

#ifdef Xyce_PARALLEL_MPI
  // Communicate singular matrix
  problem_.GetRHS()->Comm().GatherAll(&localSingularMat, &singularMat, 1);
#else
  singularMat = localSingularMat;
#endif

  if (singularMat) {

    // Inform user that singular matrix was found and linear solve has failed.
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_WARNING_0,
      "Numerically singular matrix found, returning zero solution to nonlinear solver!");

    // Put zeros in the solution since Amesos was not able to solve this problem
    problem_.GetLHS()->PutScalar( 0.0 );
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

#ifdef Xyce_VERBOSE_LINEAR
    double endSolveTime = timer_->elapsedTime();
    Xyce::lout() << "  Simple (1x1 Matrix) Solve Time: " << (endSolveTime-begSolveTime) << std::endl;
#endif

  // Output computed solution vectors, if requested.
  if (outputBaseLS_ || outputLS_) {
    if (!(base_file_number % outputBaseLS_) || !(file_number % outputLS_)) {
      char file_name[40];
      sprintf( file_name, "Simple_Soln%d.mm", base_file_number );
      EpetraExt::MultiVectorToMatrixMarketFile( file_name, *(problem_.GetLHS()) );
    }
    base_file_number++; file_number++;
  }

  // Update the total solution time
  solutionTime_ = timer_->elapsedTime();

#ifdef Xyce_VERBOSE_LINEAR
  Xyce::lout() << "Total Linear Solution Time (Simple): " << solutionTime_ << std::endl;
#endif

  return 0;
}
