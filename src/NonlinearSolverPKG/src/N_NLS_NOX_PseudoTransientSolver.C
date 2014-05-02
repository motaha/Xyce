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
// Filename       : $RCSfile: N_NLS_NOX_PseudoTransientSolver.C,v $
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
// Revision Number: $Revision: 1.19 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

/*  DEBUG: missing standard header */


#include <Xyce_config.h>


#include "N_NLS_NOX_PseudoTransientSolver.h"	// class definition
#include "N_NLS_NOX_AugmentLinSys.h"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_GlobalData.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_StatusTest_FiniteValue.H"
#include "N_NLS_LOCA_Group.h"
#include "N_ERH_ErrorMgr.h"

#include "NOX_Direction_Factory.H"
#include "NOX_LineSearch_Factory.H"
#include "NOX_Solver_SolverUtils.H"

N_NLS_NOX::PseudoTransientBased::
PseudoTransientBased(const Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys>& als,
		     const Teuchos::RefCountPtr<NOX::Abstract::Group>& xGrp, 
		     const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t, 
		     const Teuchos::RefCountPtr<Teuchos::ParameterList>& p,
		     double initialStepSize,
		     double minStepSize,
		     double maxStepSize) :
  globalData(Teuchos::rcp(new NOX::GlobalData(p))),
  augmentLSStrategy(als),
  solnPtr(xGrp),		        // pointer to xGrp
  oldSolnPtr(xGrp->clone(NOX::DeepCopy)),     // create via clone
  oldSoln(*oldSolnPtr),		        // reference to just-created pointer
  dirPtr(xGrp->getX().clone(NOX::ShapeCopy)), // create via clone 
  dir(*dirPtr),			        // reference to just-created pointer
  testPtr(t),			// pointer to t
  paramsPtr(p),		// pointer to p
  utils(*(globalData->getUtils())),                // intialize the utils
  lineSearch(NOX::LineSearch::buildLineSearch(globalData, paramsPtr->sublist("Line Search"))), // initialize the line search
  direction(NOX::Direction::buildDirection(globalData, paramsPtr->sublist("Direction"))),     // initialize the direction
  prePostOperator(globalData->getUtils(), paramsPtr->sublist("Solver Options")),
  initialStepSize_(initialStepSize),
  minStepSize_(minStepSize),
  maxStepSize_(maxStepSize),
  stepSize_(initialStepSize),
  previousStepSize_(0.0),
  group_(0),
  previousGroup_(0),
  scaleFactor_(1.0)
{
  init();
}

// Protected
void N_NLS_NOX::PseudoTransientBased::init()
{
  // Initialize 
  step_ = 0.0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;

  // Print out parameters
  if (utils.isPrintType(NOX::Utils::Parameters)) 
  {
    utils.out() << "\n" << NOX::Utils::fill(72) << "\n";
    utils.out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utils.out(),5);
  }

  checkType = NOX::Solver::parseStatusTestCheckType(paramsPtr->sublist("Solver Options"));
}

void N_NLS_NOX::PseudoTransientBased::
reset(const NOX::Abstract::Vector& initial_guess)
{
  solnPtr->setX(initial_guess);
  init();
}
 
void N_NLS_NOX::PseudoTransientBased::
reset(const NOX::Abstract::Vector& initial_guess,
      const Teuchos::RCP<NOX::StatusTest::Generic>& t)
{
  solnPtr->setX(initial_guess);
  testPtr = t;
  init();
}

N_NLS_NOX::PseudoTransientBased::~PseudoTransientBased() 
{
}


NOX::StatusTest::StatusType N_NLS_NOX::PseudoTransientBased::getStatus()
{
  return status;
}

NOX::StatusTest::StatusType N_NLS_NOX::PseudoTransientBased::step()
{
  prePostOperator.runPreIterate(*this);

  // First time thru, do some initializations
  if (this->getNumIterations() == 0) {
    
    // Compute F of initital guess
    NOX::Abstract::Group::ReturnType rtype = solnPtr->computeF();
    if (rtype != NOX::Abstract::Group::Ok) {
      Xyce::dout() << "NOX::Solver::PseudoTransientBased::step - Unable to compute F" << std::endl;
      throw "NOX Error";
    }
    
    // Test the initial guess
    status = testPtr->checkStatus(*this, checkType);
    if ((status == NOX::StatusTest::Converged) &&
	(utils.isPrintType(NOX::Utils::Warning))) {
      utils.out() << "Warning: NOX::Solver::PseudoTransientBased::step() - The solution passed "
		  << "into the solver (either through constructor or reset method) "
		  << "is already converged!  The solver will not "
		  << "attempt to solve this system since status is flagged as "
		  << "converged." << std::endl;
    }
    
    // Print out status tests
    if (utils.isPrintType(NOX::Utils::Parameters)) {
      utils.out() << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
      testPtr->print(utils.out(), 5);
      utils.out() <<"\n" << NOX::Utils::fill(72) << "\n";
    }
    
    if (status != NOX::StatusTest::Unconverged) {
      prePostOperator.runPostIterate(*this);
      return status;
    }
  }

  // Compute a new stp length
  if (this->getNumIterations() == 0) {
    stepSize_ = initialStepSize_;
  }
  else {
    previousStepSize_ = stepSize_;
    double f_n = this->getSolutionGroup().getNormF();
    double f_nm1 = this->getPreviousSolutionGroup().getNormF();
    
    stepSize_ = scaleFactor_ * previousStepSize_ * f_nm1 / f_n;

    //Xyce::dout() << "f_nm1/f_n = " << f_nm1 / f_n << std::endl;

    //if ((stepSize_ > previousStepSize_) && (f_n > 1.0e+10)) {
    //stepSize_ = previousStepSize_;
    //}
    
    //if (stepSize_ < previousStepSize_ * 0.1 )
    //stepSize_ = previousStepSize_ * 0.1;
    
    //if (stepSize_ > previousStepSize_ * 1.0e4 )
    //stepSize_ = previousStepSize_ * 10.0;

    if (stepSize_ < minStepSize_)
      stepSize_ = minStepSize_;
    
    if (stepSize_ > maxStepSize_)
      stepSize_ = maxStepSize_;
    
  }

  // Set Steplength with group
//   group_->setPseudoTransientStepSize(stepSize_);
//   previousGroup_->setPseudoTransientStepSize(stepSize_);
  augmentLSStrategy->setProgressVariable(stepSize_);
  
  // First check status
  if (status != NOX::StatusTest::Unconverged) {
    prePostOperator.runPostIterate(*this);
    return status;
  }

  // Copy pointers into temporary references
  NOX::Abstract::Group& soln = *solnPtr;
  NOX::StatusTest::Generic& test = *testPtr;

  // Compute the direction for the update vector at the current solution.
  bool ok;
  ok = direction->compute(dir, soln, *this);
  if (!ok) 
  {
    Xyce::dout() << "N_NLS_NOX::PseudoTransientBased::iterate - unable to calculate direction" << std::endl;
    status = NOX::StatusTest::Failed;
    prePostOperator.runPostIterate(*this);
    return status;
  }

  // Update iteration count.
  nIter ++;

  // Copy current soln to the old soln.
  oldSoln = soln;

  // Do line search and compute new soln.
  ok = lineSearch->compute(soln, step_, dir, *this);
  if (!ok) 
  {
    if (nIter == 0) 
    {
      Xyce::dout() << "N_NLS_NOX::PseudoTransientBased::iterate - line search failed" << std::endl;
      status = NOX::StatusTest::Failed;
      prePostOperator.runPostIterate(*this);
      return status;
    }
    else if (utils.isPrintType(NOX::Utils::Warning))
      utils.out() << "N_NLS_NOX::PseudoTransientBased::iterate - using recovery step for line search" << std::endl;
  }

  // Compute F for new current solution.
  NOX::Abstract::Group::ReturnType rtype = soln.computeF();
  if (rtype != NOX::Abstract::Group::Ok) 
  {
    utils.out() << "N_NLS_NOX::PseudoTransientBased::iterate - unable to compute F" << std::endl;
    status = NOX::StatusTest::Failed;
    prePostOperator.runPostIterate(*this);
    return status;
  }
  
  NOX::StatusTest::FiniteValue fv;
  NOX::StatusTest::StatusType fvStatus = fv.checkStatus(*this, checkType);
  if (fvStatus == NOX::StatusTest::Failed) {
    scaleFactor_ *= 0.5;
    *group_ = *previousGroup_;
    prePostOperator.runPostIterate(*this);
    group_->setX(group_->getX());
    if (stepSize_ > minStepSize_)
      return (NOX::StatusTest::Unconverged);
    else
      return (NOX::StatusTest::Failed);
  }
  else 
    scaleFactor_ = 1.0;

  // Evaluate the current status.
  status = test.checkStatus(*this, checkType);
 
  prePostOperator.runPostIterate(*this);

  // Return status.
  return status;
}

NOX::StatusTest::StatusType N_NLS_NOX::PseudoTransientBased::solve()
{
  prePostOperator.runPreSolve(*this);

  
  group_ = 0;
  NOX::Abstract::Group* tmpGroup = 
    const_cast<NOX::Abstract::Group*>(&(this->getSolutionGroup()));
  group_ = dynamic_cast<N_NLS_LOCA::Group*>(tmpGroup);

  
  previousGroup_ = 0;
  NOX::Abstract::Group* tmpPreviousGroup = 
    const_cast<NOX::Abstract::Group*>(&(this->getPreviousSolutionGroup()));
  previousGroup_ = dynamic_cast<N_NLS_LOCA::Group*>(tmpPreviousGroup);

  if ((group_ == 0) || (previousGroup_ == 0)) {
    std::string message = "N_NLS_NOX::PrePostOperator\n Failed to dynamic_cast the old and new groups to N_NLS_LOCA::Groups! ";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

  // Enable xyce group pseudotransient
  group_->setAugmentLinearSystem(true, augmentLSStrategy);
  previousGroup_->setAugmentLinearSystem(true, augmentLSStrategy);
    
  printUpdate();

  // Iterate until converged or failed
  while (status == NOX::StatusTest::Unconverged) 
  {
    status = step();
    printUpdate();
  }

  group_->setAugmentLinearSystem(false, augmentLSStrategy);
  previousGroup_->setAugmentLinearSystem(false, augmentLSStrategy);

  Teuchos::ParameterList& outputParams = paramsPtr->sublist("Output");
  outputParams.set("Nonlinear Iterations", nIter);
  outputParams.set("2-Norm of Residual", solnPtr->getNormF());

  prePostOperator.runPostSolve(*this);

  return status;
}

const NOX::Abstract::Group& 
N_NLS_NOX::PseudoTransientBased::getSolutionGroup() const
{
  return *solnPtr;
}

const NOX::Abstract::Group& 
N_NLS_NOX::PseudoTransientBased::getPreviousSolutionGroup() const
{
  return oldSoln;
}

int N_NLS_NOX::PseudoTransientBased::getNumIterations() const
{
  return nIter;
}

const Teuchos::ParameterList& 
N_NLS_NOX::PseudoTransientBased::getList() const
{
  return *paramsPtr;
}

// protected
void N_NLS_NOX::PseudoTransientBased::printUpdate() 
{
  double normSoln = 0;
  double normStep = 0;

  // Print the status test parameters at each iteration if requested  
  if ((status == NOX::StatusTest::Unconverged) && 
      (utils.isPrintType(NOX::Utils::OuterIterationStatusTest))) 
  {
    utils.out() << NOX::Utils::fill(72) << "\n";
    utils.out() << "-- Status Test Results --\n";    
    testPtr->print(utils.out());
    utils.out() << NOX::Utils::fill(72) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (utils.isPrintType(NOX::Utils::OuterIteration)) 
  {
    normSoln = solnPtr->getNormF();
    normStep = (nIter > 0) ? dir.norm() : 0;
  }

  // ...But only the print process actually prints the result.
  if (utils.isPrintType(NOX::Utils::OuterIteration)) 
  {
    utils.out() << "\n" << NOX::Utils::fill(72) << "\n";
    utils.out() << "-- Nonlinear Solver Step " << nIter << " -- \n";
    utils.out() << "f = " << utils.sciformat(normSoln);
    utils.out() << "  step = " << utils.sciformat(step_);
    utils.out() << "  dx = " << utils.sciformat(normStep);
    if (status == NOX::StatusTest::Converged)
      utils.out() << " (Converged!)";
    if (status == NOX::StatusTest::Failed)
      utils.out() << " (Failed!)";
    utils.out() << "\n" << NOX::Utils::fill(72) << "\n" << std::endl;
  }

  // Print the final parameter values of the status test
  if ((status != NOX::StatusTest::Unconverged) && 
      (utils.isPrintType(NOX::Utils::OuterIteration))) 
  {
    utils.out() << NOX::Utils::fill(72) << "\n";
    utils.out() << "-- Final Status Test Results --\n";    
    testPtr->print(utils.out());
    utils.out() << NOX::Utils::fill(72) << "\n";
  }
}

double N_NLS_NOX::PseudoTransientBased::getStepSize() const
{
  return step_;
}

double N_NLS_NOX::PseudoTransientBased::getPseudoTransientStepSize() const
{
  return stepSize_;
}

