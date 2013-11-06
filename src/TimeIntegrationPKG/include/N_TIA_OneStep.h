//-----------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: N_TIA_OneStep.h,v $
//
// Purpose       : This file defines the classes for trapezoidal method 
//                 order 1-2 method.
//
// Special Notes :
//
// Creator       : Ting Mei
//
// Creation Date : 10/31/07
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.14.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_ONE_STEP_H
#define Xyce_N_TIA_ONE_STEP_H


// ---------- Standard Includes ----------
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_PDS_ParMap.h>

//-----------------------------------------------------------------------------
// Class         : N_TIA_OneStep
// Purpose       : OneStep Formula Integration Class
//                 (derived from N_TIA_TimeIntegrationMethod)
// Special Notes :
// Creator       : Ting Mei, SNL
// Creation Date : 10/31/07
//-----------------------------------------------------------------------------
class N_TIA_OneStep : public N_TIA_TimeIntegrationMethod
{
public:

  // Destructor
  ~N_TIA_OneStep() {}

  // Predict solution at next time point.
  virtual void obtainPredictor();

  // Evaluate the predictor derivative formula.
  virtual void obtainPredictorDeriv() {return;}

  // Evaluate the corrector derivative formula.
  virtual void obtainCorrectorDeriv() {return;}

  // Compute an estimate of the error in the integration step.
  virtual void updateDerivsBlock(const list < index_pair > & solGIDList,
                                 const list < index_pair > & staGIDList) {return;}

  // Compute an estimate of the error in the integration step.
  virtual double computeErrorEstimate() { return sec.ck_*ds.WRMS_errorNorm(); }

  // Interpolate solution approximation at prescribed time point.
  virtual bool interpolateSolution(double timepoint,
                  N_LAS_Vector * tmpSolVectorPtr, vector<N_LAS_Vector*> & historyVec);

  // Interpolate MPDE solution approximation at prescribed time point.
  virtual bool interpolateMPDESolution(vector<double>& timepoint,
                  N_LAS_Vector * tmpSolVectorPtr);
  
  // Print transient output from MPDE simulation
  virtual bool printMPDEOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  const double time,
	                N_LAS_Vector * solnVecPtr,
                  const std::vector<double> & fastTimes );

  // Print transient output from WaMPDE simulation
  virtual bool printWaMPDEOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  const double time,
	                N_LAS_Vector * solnVecPtr,
                  const std::vector<double> & fastTimes,
                  const int phiGID );

  // Print output using interpolation when order is high
  virtual bool printOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  const double time,
	                N_LAS_Vector * solnVecPtr,
                  const bool doNotInterpolate,
                  const vector<double> & outputInterpolationTimes,
                  bool skipPrintLineOutput );

  // .SAVE output using interpolation when order is high
  virtual bool saveOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  N_LAS_Vector * solnVecPtr,
                  const double saveTime,
                  const bool doNotInterpolate); 

  // Computes the step adjustment.
  virtual double computeExpoStepAdjust(double stepadjust);

  // Gets the time-integration order.
  virtual int getOrder() {return sec.currentOrder_;}
  virtual int getUsedOrder() {return sec.usedOrder_;}
  virtual int getNumberOfSteps () { return sec.numberOfSteps_; }
  virtual int getNscsco () { return sec.nscsco_; }

  virtual void getInitialQnorm (N_TIA_TwoLevelError & tle);
  virtual void setupTwoLevelError(N_TIA_TwoLevelError & tle);

  // 2/16/04 tscoffe:  This is just for testing, we'll need a dynamic fcn later

  // Time-integration factory.
  static N_TIA_TimeIntegrationMethod * factory
    (N_TIA_TIAParams & tiaP, 
     N_TIA_StepErrorControl & secTmp ,
     N_TIA_DataStore & dsTmp);

  // 03/08/04 tscoffe:  New functions necessary for OneStep
  // Evaluate corrector residual for nonlinear solver
  virtual void obtainResidual();

  // Evaluate corrector Jacobian for nonlinear solver
  virtual void obtainJacobian();

  // Update history array after a successful step 
  virtual void updateHistory();

  // Restore history array after a failed step
  virtual void restoreHistory();

  // Return max order of method (this should obey user option maxorder)
  virtual int getMaxOrder(){ return sec.maxOrder_;}

  // Update method coefficients
  virtual void updateCoeffs();

  // Initialize method with initial solution & step-size
  virtual void initialize();

  // setup 2-level data.
  virtual void setTwoLevelTimeInfo(const N_TIA_TimeIntInfo & tiInfo);

  // Reject a step (this turns back history and chooses new order & step-size)
  virtual void rejectStep();

  // Reject a step, but only restore the history, as the new step will be
  // imposed by an outside program. 
  virtual void rejectStepForHabanero ();

  // Complete a step (this updates history and chooses new order & step-size)
  virtual void completeStep();

  virtual void updateStateDeriv ();
  
  // calculates dQ/dt component of store vector and adds it to store vector
  virtual void updateLeadCurrent ();

private:

  // Constructor (private).
  N_TIA_OneStep(N_TIA_TIAParams & tiaP,
                                  N_TIA_StepErrorControl & secTmp,
                                  N_TIA_DataStore & dsTmp);

  // Check whether to reduce order independent of local error test
  virtual void checkReduceOrder();

  enum actionFlag { TIAAction_UNSET, TIAAction_LOWER, TIAAction_MAINTAIN, TIAAction_RAISE };

  // Used to keep track of last interpolation point in printMPDEOutputSolution
  double timept_; 
  
  // used for second order interpolation where the last two times steps are needed.
  // lastTimeStep is stored in StepErrorControl and the last, last time step is stored
  // here as it's only needed by trap.  Called timeStepForHistory2 as it's the 
  // time step associated with the step in xHistory[2]
  double timeStepForHistory2_;
};

  // Computes the step addjustment.
  // 2/16/04 tscoffe:  I'm not exactly sure what this routine is for...
inline double N_TIA_OneStep::computeExpoStepAdjust(
  double stepadjust)
{
  return pow(stepadjust, 1.0 / 3.0);
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_OneStep::factory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :  10/31/07 
//-----------------------------------------------------------------------------
inline N_TIA_TimeIntegrationMethod * N_TIA_OneStep::factory
  (N_TIA_TIAParams & tiaP, 
   N_TIA_StepErrorControl & secTmp ,
   N_TIA_DataStore & dsTmp)
{
  N_TIA_TimeIntegrationMethod * timPtr =
    new N_TIA_OneStep(tiaP,secTmp,dsTmp);
  return timPtr;
}

#endif     //Xyce_N_TIA_ONE_STEP_H

