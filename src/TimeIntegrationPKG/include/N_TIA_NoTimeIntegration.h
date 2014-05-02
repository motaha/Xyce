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

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: N_TIA_NoTimeIntegration.h,v $
//
// Purpose       : This file defines the classes for the "no time integration"
//                 method.
//
// Special Notes :
//
// Creator       : Eric Keiter
//
// Creation Date : 7/21/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.32 $
//
// Revision Date  : $Date: 2014/02/24 23:49:26 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_NO_TIME_INTEGRATION_H
#define Xyce_N_TIA_NO_TIME_INTEGRATION_H

// ---------- Standard Includes ----------
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TimeIntegrationMethods.h>

#define MATRIX_FAILSAFE 1
#define NUM_LIMIT  1.0e-20

//-----------------------------------------------------------------------------
// Class         : N_TIA_NoTimeIntegration
// Purpose       : Class objects for use during DC analysis (applies only when
//                 all time derivatives are set to 0) (derived from
//                 N_TIA_TimeIntegrationMethod)
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class N_TIA_NoTimeIntegration : public N_TIA_TimeIntegrationMethod
{
public:

  // Destructor
  ~N_TIA_NoTimeIntegration();

  // Predict solution at next time point (No Integration).
  virtual void obtainPredictor() { ds.usePreviousSolAsPredictor (); }

  // Evaluate the predictor derivative formula (No Integration).
  virtual void obtainPredictorDeriv() { }

  // Evaluate the corrector derivative formula (No Integration).
  virtual void obtainCorrectorDeriv();

  // Compute an estimate of the error in the integration step (No Integration).
  virtual void updateDerivsBlock(const std::list< index_pair > & solGIDList,
                                 const std::list< index_pair > & staGIDList)
    { }

#ifdef MATRIX_FAILSAFE
  virtual double partialTimeDeriv() { return NUM_LIMIT; }

#else
  virtual double partialTimeDeriv() { return 0.0; }

#endif

  // Compute an estimate of the error in the integration step (No Integration).
  virtual double computeErrorEstimate() { return 0.0; }

  // Interpolate solution approximation at prescribed time point (No
  // Integration).
  virtual bool interpolateSolution(double timepoint,
      	N_LAS_Vector * tmpSolVectorPtr, std::vector<N_LAS_Vector*> & historyVec ) 
    {return false;};

  // Computes the step adjustment (No Integration).
  virtual double computeExpoStepAdjust(double stepadjust) { return 0.0; }

  // Gets the time-integration order (No Integration).
  virtual int getOrder() { return 0; }
  virtual int getUsedOrder() { return 0; }

  virtual void getInitialQnorm (N_TIA_TwoLevelError & tle);
  virtual void setupTwoLevelError(N_TIA_TwoLevelError & tle);

  // Time-integration factory.
  static N_TIA_TimeIntegrationMethod * factory(N_TIA_TIAParams & tiaP, 
                                               N_TIA_StepErrorControl & secTmp ,
                                               N_TIA_DataStore & dsTmp);

  // 03/08/04 erkeite:  New functions necessary new-DAE:
  // Evaluate residual for nonlinear solver
  void obtainResidual();

  // Evaluate Jacobian for nonlinear solver
  void obtainJacobian();

  // Apply Jacobian for nonlinear solver
  void applyJacobian(const N_LAS_Vector& input, N_LAS_Vector& result);

private:
  // Copied from N_TIA_BackwardDifferentiation15.h
  double alphas;    // $\alpha_s$ fixed-leading coefficient of this BDF method

  // Constructor (private).
  N_TIA_NoTimeIntegration(N_TIA_TIAParams & tiaP,
                          N_TIA_StepErrorControl & secTmp,
                          N_TIA_DataStore & dsTmp); 
};

//-----------------------------------------------------------------------------
// Function      : N_TIA_NoTimeIntegration::factory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/08/01
//-----------------------------------------------------------------------------
inline N_TIA_TimeIntegrationMethod * N_TIA_NoTimeIntegration::factory
   (N_TIA_TIAParams & tiaP, 
    N_TIA_StepErrorControl & secTmp ,
    N_TIA_DataStore & dsTmp) 
{
  N_TIA_TimeIntegrationMethod * timPtr = 
    new N_TIA_NoTimeIntegration(tiaP,secTmp,dsTmp);
  return timPtr;
}

#endif // Xyce_N_TIA_NO_TIME_INTEGRATION_H
