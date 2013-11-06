//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_LOCA_StepSizeControl.C,v $
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
// Revision Number: $Revision: 1.3.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:48 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

/* DEBUG:  missing standard header */


#include <Xyce_config.h>



#include "N_NLS_LOCA_StepSizeControl.h"

#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_Stepper.H"
#include "NOX_Solver_Generic.H"
#include "LOCA_Utils.H"
#include "LOCA_MultiContinuation_AbstractStrategy.H"
#include "LOCA_MultiContinuation_ExtendedVector.H"
#include "LOCA_NewStepper.H"

N_NLS_LOCA::StepSizeControl::StepSizeControl() :
  maxStepSize(0.0),
  minStepSize(0.0),
  agrValue(0.0)
{

}

N_NLS_LOCA::StepSizeControl::~StepSizeControl()
{
}

NOX::Abstract::Group::ReturnType 
N_NLS_LOCA::StepSizeControl::reset(NOX::Parameter::List& params) 
{
  maxStepSize = params.getParameter("Max Step Size", 1.0e+12);
  minStepSize = params.getParameter("Min Step Size", 1.0e-12);
  startStepSize = params.getParameter("Initial Step Size", 1.0);
  failedFactor = params.getParameter("Failed Step Reduction Factor", 0.5);
  successFactor = params.getParameter("Successful Step Increase Factor", 1.26);
  prevStepSize = 0.0;
  isFirstStep = true;  
  agrValue = params.getParameter("Aggressiveness", 0.0);

  return NOX::Abstract::Group::Ok;
} 

NOX::Abstract::Group::ReturnType 
N_NLS_LOCA::StepSizeControl::compute(
		       LOCA::Continuation::ExtendedGroup& curGroup,
		       const LOCA::Continuation::ExtendedVector& predictor,
		       const NOX::Solver::Generic& solver,
		       const LOCA::Abstract::Iterator::StepStatus& stepStatus,
		       const LOCA::Stepper& stepper,
		       double& stepSize) 
{
  // If this is the first step, set step size to initial value
  if (isFirstStep) {
    double dpds = predictor.getParam();
    if (dpds != 0.0) {
      startStepSize /= dpds;
      maxStepSize /= dpds;
      minStepSize /= dpds;
    }
    isFirstStep = false;
    stepSize = startStepSize;
    prevStepSize = 0.0;
  }
  else {
    
    // A failed nonlinear solve cuts the step size in half
    if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful) {
      stepSize *= failedFactor;    
    }
    else {

      double ds_ratio = curGroup.getStepSizeScaleFactor();
      startStepSize *= ds_ratio;
      maxStepSize *= ds_ratio;
      minStepSize *= ds_ratio;

      // Get maximum number of nonlinear iterations from stepper parameters
      const NOX::Parameter::List& p = LOCA::Utils::getSublist("Stepper");
      double maxNonlinearSteps 
	= static_cast<double>(p.getParameter("Max Nonlinear Iterations", 15));
      
      // Get number of nonlinear iterations in last step
      double numNonlinearSteps = 
	static_cast<double>(solver.getNumIterations());

      // Save successful stepsize as previous
      prevStepSize = stepSize;

      // adapive step size control
      double factor = (maxNonlinearSteps - numNonlinearSteps) 
               	      / (maxNonlinearSteps);

      stepSize *= (1.0 + agrValue * factor * factor);

      stepSize *= ds_ratio;
    } 
  }

  // Clip step size to be within prescribed bounds
  NOX::Abstract::Group::ReturnType res = clipStepSize(stepSize);

  return res;
}

NOX::Abstract::Group::ReturnType 
N_NLS_LOCA::StepSizeControl::compute(
		     LOCA::MultiContinuation::AbstractStrategy& curGroup,
		     const LOCA::MultiContinuation::ExtendedVector& predictor,
		     const NOX::Solver::Generic& solver,
		     const LOCA::Abstract::Iterator::StepStatus& stepStatus,
		     const LOCA::NewStepper& stepper,
		     double& stepSize) 
{
  // If this is the first step, set step size to initial value
  if (isFirstStep) {
    double dpds = predictor.getScalar(0);
    if (dpds != 0.0) {
      startStepSize /= dpds;
      maxStepSize /= dpds;
      minStepSize /= dpds;
    }
    isFirstStep = false;
    stepSize = startStepSize;
    prevStepSize = 0.0;
  }
  else {
  
    // A failed nonlinear solve cuts the step size in half
    if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful) {
      stepSize *= failedFactor;    
    }
    else {

      double ds_ratio = curGroup.getStepSizeScaleFactor();
      startStepSize *= ds_ratio;
      maxStepSize *= ds_ratio;
      minStepSize *= ds_ratio;

      // Get maximum number of nonlinear iterations from stepper parameters
      const NOX::Parameter::List& p = LOCA::Utils::getSublist("Stepper");
      double maxNonlinearSteps 
	= static_cast<double>(p.getParameter("Max Nonlinear Iterations", 15));
      
      // Get number of nonlinear iterations in last step
      double numNonlinearSteps = 
	static_cast<double>(solver.getNumIterations());

      // Save successful stepsize as previous
      prevStepSize = stepSize;

      // adapive step size control
      double factor = (maxNonlinearSteps - numNonlinearSteps) 
               	      / (maxNonlinearSteps);

      stepSize *= (1.0 + agrValue * factor * factor);

      stepSize *= ds_ratio;
    } 
  }

  // Clip step size to be within prescribed bounds
  NOX::Abstract::Group::ReturnType res = clipStepSize(stepSize);

  return res;
}


NOX::Abstract::Group::ReturnType
N_NLS_LOCA::StepSizeControl::clipStepSize(double& stepSize)
{
  NOX::Abstract::Group::ReturnType res = NOX::Abstract::Group::Ok;

  // Compute sign of step size
  double signStep = 1.0;
  if (stepSize < 0.0)
    signStep = -1.0;

  // Clip the step size if above the bounds
  if (fabs(stepSize) > maxStepSize)
    stepSize = signStep*maxStepSize;

  // Clip step size at minimum, signal for failed run
  if (fabs(stepSize) < minStepSize) {
    res = NOX::Abstract::Group::Failed;
    stepSize =  signStep*minStepSize;
    if (LOCA::Utils::doPrint(LOCA::Utils::Error)) {
      cout << "\n\tStep size reached minimum step size bound"
           << endl;
    }
  }

  return res;
}

double
N_NLS_LOCA::StepSizeControl::getPrevStepSize() const {
  return prevStepSize;
}

double
N_NLS_LOCA::StepSizeControl::getStartStepSize() const {
  return startStepSize;
}

