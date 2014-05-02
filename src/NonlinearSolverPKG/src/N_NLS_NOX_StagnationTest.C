//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_NOX_StagnationTest.C,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : 
//
// Creation Date  : 
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

/* DEBUG:  missing standard header */


#include <Xyce_config.h>


#include "N_NLS_NOX_StagnationTest.h" // class definition
#include "NOX_Common.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Abstract_Group.H"

N_NLS_NOX::Stagnation::Stagnation(int maxSteps_, double tolerance_) :
  maxSteps(maxSteps_),
  numSteps(0),
  lastIteration(-1),
  tolerance(tolerance_),
  convRate(1.0),
  minConvRate(1.0),
  normFInit(1.0),
  status(NOX::StatusTest::Unconverged)
{
    
}

N_NLS_NOX::Stagnation::~Stagnation()
{

}

NOX::StatusTest::StatusType 
N_NLS_NOX::Stagnation::checkStatus(const NOX::Solver::Generic& problem)
{
  status = NOX::StatusTest::Unconverged;

  // First time through we don't do anything but reset the counters
  int niters = problem.getNumIterations(); 
  if (niters == 0) {
    numSteps = 0;
    lastIteration = 0;
    convRate = 1.0;
    //minConvRate = 1.0;  // Don't reset this.  Xyce solver never does.
    normFInit = problem.getSolutionGroup().getNormF();
    return NOX::StatusTest::Unconverged;
  } 

  // Make sure we have not already counted the last nonlinear iteration.
  // This protects against multiple calls to checkStatus() in between 
  // nonlinear iterations.
  bool isCounted = false;
  if (niters == lastIteration) {
    isCounted = true;
  }
  else
    lastIteration = niters;

  // Compute the convergenc rate and set counter appropriately
  if (!isCounted) {

    convRate = problem.getSolutionGroup().getNormF() / 
               problem.getPreviousSolutionGroup().getNormF();
    
    if (fabs(convRate - 1.0) <= tolerance) {
      
      if ((numSteps == 0) || (convRate < minConvRate)) 
	minConvRate = convRate;
      
      ++numSteps ;
    }
    else
      numSteps = 0;
   
  }

  if (numSteps >= maxSteps) {
  
    double initConvRate = problem.getSolutionGroup().getNormF()/normFInit;

    if ((initConvRate <= 0.9) && (minConvRate <= 1.0)) {
      status = NOX::StatusTest::Converged;
    }
    else
      status = NOX::StatusTest::Failed;

  }

  return status;
}

NOX::StatusTest::StatusType N_NLS_NOX::Stagnation::getStatus() const
{
  return status;
}

std::ostream& N_NLS_NOX::Stagnation::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << status;
  stream << "Stagnation Count = " << numSteps << " = " << maxSteps << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "             (convergence rate = " << convRate << ")";
  
  if (status == NOX::StatusTest::Converged) {
    stream << "\n             (Near Convergence!)" << std::endl;
  }
  
  stream << std::endl;
 return stream;
}


int N_NLS_NOX::Stagnation::getMaxNumSteps() const
{
  return maxSteps;
}

int N_NLS_NOX::Stagnation::getCurrentNumSteps() const
{
  return numSteps;
}

double N_NLS_NOX::Stagnation::getTolerance() const
{
  return tolerance;
}

double N_NLS_NOX::Stagnation::getConvRate() const
{
  return convRate;
}

