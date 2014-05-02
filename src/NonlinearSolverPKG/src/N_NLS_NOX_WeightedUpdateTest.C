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
// Filename       : $RCSfile: N_NLS_NOX_WeightedUpdateTest.C,v $
//
// Purpose        : Status test.
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.17 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include "N_NLS_NOX_WeightedUpdateTest.h"
#include "N_NLS_NOX_Vector.h"
#include "N_LAS_Vector.h"
#include "NOX.H"
#include "NOX_Solver_LineSearchBased.H"

// ----------   Namespaces   ----------

using namespace N_NLS_NOX;

// ----------   Code   ----------

WeightedUpdateTest::WeightedUpdateTest(N_LAS_Vector** currSolVectorPtrPtr,
				       double epsilon_a, double epsilon_r, 
				       double tol, bool isTransient) :
  epsilon_a_(epsilon_a),
  epsilon_r_(epsilon_r),
  tol_(tol),
  isTransient_(isTransient)
{
  result_ = 0.0;
  status_ = NOX::StatusTest::Unconverged;
  oldTimeStepVectorPtrPtr_ = currSolVectorPtrPtr;
  weightsVectorPtr_ = 0;
  updateVectorPtr_ = 0;
  tmpVectorPtr_ = 0;
}

WeightedUpdateTest::~WeightedUpdateTest()
{
  if( weightsVectorPtr_ != 0 )
  {
    delete weightsVectorPtr_;
    delete updateVectorPtr_;
    delete tmpVectorPtr_;
  }
}


NOX::StatusTest::StatusType WeightedUpdateTest::checkStatus(const NOX::Solver::Generic& problem)
{
  status_ = NOX::StatusTest::Unconverged;

  int niters = problem.getNumIterations();

  // Copy into local reference
  N_LAS_Vector& oldTimeStepX = **oldTimeStepVectorPtrPtr_;

  const N_LAS_Vector& x = dynamic_cast<const N_NLS_NOX::Vector&>
    (problem.getSolutionGroup().getX()).getNativeVectorRef();
  const N_LAS_Vector& oldX = dynamic_cast<const N_NLS_NOX::Vector&>
    (problem.getPreviousSolutionGroup().getX()).getNativeVectorRef();

  // Allocate space if necessary
  if (weightsVectorPtr_ == 0) {
    weightsVectorPtr_ = new N_LAS_Vector(x);
    updateVectorPtr_ = new N_LAS_Vector(x);
    tmpVectorPtr_ = new N_LAS_Vector(x);
  }

  // Local references
  N_LAS_Vector& weights = *weightsVectorPtr_;
  N_LAS_Vector& update = *updateVectorPtr_;
  N_LAS_Vector& tmp = *tmpVectorPtr_;

  // Compute local portion of weights vector
  // Weights are recomputed at each nonlinear iteration of a DC Op calc
  // but only at the beginning of a transient nonlinear solve.
  if ((!isTransient_) || (niters == 0)) {
    int length = x.localLength();
    for (int i = 0; i < length; ++i ) {
      //update[i] = x[i] - oldX[i];
      weights[i] =
        epsilon_r_ * Xycemax(fabs(x[i]), fabs(oldTimeStepX[i])) +
        epsilon_a_;
    }
    //weights.printPetraObject();
  }

  if (niters < 1) {
    result_ = 0.0;
    return status_;
  }

  // Next compute the update
  update.update(1.0, x, -1.0, oldX, 0.0);

  // Compute final result
#ifdef Xyce_SPICE_NORMS
  update.wMaxNorm(weights,tmp,&result_);
#else
  update.wRMSNorm(weights,&result_);
#endif

  //cout.precision(12);
  //cout << "\n ** WRMS = "<<result_<< "\n" <<endl;
 
  // RPP: If a line search is being used, we must account for any 
  // damping of the step length.  Otherwise delta X could be small due 
  // the line search and not due to being close to a solution.
  const NOX::Solver::LineSearchBased* test = 0;
  test = dynamic_cast<const NOX::Solver::LineSearchBased*>(&problem);
  if (test != 0) {
    result_ = result_/(test->getStepSize());
  }

  if (result_ < tol_)
    status_ = NOX::StatusTest::Converged;

  return status_;
}

std::ostream& WeightedUpdateTest::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << status_;
  stream << "Weighted Update = " << NOX::Utils::sciformat(result_, 3) << " < " << NOX::Utils::sciformat(tol_, 3) << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << std::setw(13) << " ";
  stream << "(with e_r = " << NOX::Utils::sciformat(epsilon_r_, 3);
  stream << " and e_a = " << NOX::Utils::sciformat(epsilon_a_, 3) << ")"; ;
  stream << std::endl;
  return stream;
}

