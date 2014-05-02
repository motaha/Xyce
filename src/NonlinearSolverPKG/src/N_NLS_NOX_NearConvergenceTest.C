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
// Filename       : $RCSfile: N_NLS_NOX_NearConvergenceTest.C,v $
//
// Purpose        : Status test.
//
// Special Notes  :
//
// Creator        : Roger Pawlowski, SNL 9233
//
// Creation Date  : 04/02/2003
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include "N_NLS_NOX_NearConvergenceTest.h"
#include "N_NLS_NOX_Vector.h"
#include "N_LAS_Vector.h"
#include "NOX.H"
#include "NOX_Solver_LineSearchBased.H"

// ----------   Namespaces   ----------

using namespace N_NLS_NOX;

// ----------   Code   ----------

NearConvergenceTest::NearConvergenceTest(int maxIters,
					 double convRate,
					 double relativeConvRate) :
  maxIters_(maxIters),
  nIters_(0),
  requestedConvRate_(convRate),
  currentConvRate_(1.0),
  requestedRelativeConvRate_(relativeConvRate),
  currentRelativeConvRate_(1.0),
  normResidualInit_(1.0),
  status_(NOX::StatusTest::Unconverged)
{

}

NearConvergenceTest::~NearConvergenceTest()
{

}


NOX::StatusTest::StatusType NearConvergenceTest::checkStatus(const NOX::Solver::Generic& problem)
{
  status_ = NOX::StatusTest::Unconverged;

  nIters_ = problem.getNumIterations();

  if (nIters_ == 0)
    normResidualInit_ = problem.getSolutionGroup().getNormF();
  
  // Test only if we hit the max number of iterations
  if (nIters_ >= maxIters_) {
    
    // ||F(x_current)|| / ||F(x_previous)||
    currentConvRate_ = (problem.getSolutionGroup().getNormF()) / 
                       (problem.getPreviousSolutionGroup().getNormF());

    //  ||F(x)|| / ||F(x_init)||
    currentRelativeConvRate_ = (problem.getSolutionGroup().getNormF()) / 
                                   normResidualInit_;
    
    if ((currentConvRate_ <= requestedConvRate_) && 
	(currentRelativeConvRate_ <= requestedRelativeConvRate_)) {
      status_ = NOX::StatusTest::Converged;
    }
    else
      status_ = NOX::StatusTest::Failed;
  }

  return status_;
}

std::ostream& NearConvergenceTest::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; ++j )
    stream << ' ';

  stream << status_;

  stream << "Near Convergence = ";

  if (status_ == NOX::StatusTest::Converged)
    stream << "Yes\n";
  else
    stream << "No\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "             (Max Iters = " << nIters_ << " = " << maxIters_ << ")\n";  
  
  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "             (Conv Rate = " << NOX::Utils::sciformat(currentConvRate_,3) << " < " << NOX::Utils::sciformat(requestedConvRate_,3) << ")\n";  
    
  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "             (Rel Conv Rate = " << NOX::Utils::sciformat(currentRelativeConvRate_, 3) << " < " << NOX::Utils::sciformat(requestedRelativeConvRate_, 3) << ")";  

  stream << std::endl;

  return stream;
}

