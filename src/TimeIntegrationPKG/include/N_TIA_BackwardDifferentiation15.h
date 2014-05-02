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
// Filename      : $RCSfile: N_TIA_BackwardDifferentiation15.h,v $
//
// Purpose       : This file defines the classes for the backward
//                 differentiation, order 1-5 method.
//
// Special Notes :
//
// Creator       : Todd Coffey
//
// Creation Date : 2/16/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.51 $
//
// Revision Date  : $Date: 2014/02/24 23:49:26 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_BACKWARD_DIFFERENTIATION_15_H
#define Xyce_N_TIA_BACKWARD_DIFFERENTIATION_15_H


// ---------- Standard Includes ----------
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

#include <N_ANP_fwd.h>

#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_PDS_ParMap.h>

//-----------------------------------------------------------------------------
// Class         : N_TIA_BackwardDifferentiation15
// Purpose       : Backward Differentiation Formula (of order 1-5)
//                 Integration Class
//                 (derived from N_TIA_TimeIntegrationMethod)
// Special Notes :
// Creator       : Todd Coffey, SNL
// Creation Date : 2/16/04
//-----------------------------------------------------------------------------
class N_TIA_BackwardDifferentiation15 : public N_TIA_TimeIntegrationMethod
{
public:

  // Destructor
  ~N_TIA_BackwardDifferentiation15() {}

  // Predict solution at next time point (BDF1-5).
  virtual void obtainPredictor();

  // Evaluate the predictor derivative formula (BDF1-5).
  virtual void obtainPredictorDeriv() {return;}

  // Evaluate the corrector derivative formula (BDF1-5).
  virtual void obtainCorrectorDeriv() {return;}

  // Compute an estimate of the error in the integration step (BDF1-5).
  virtual void updateDerivsBlock(const std::list< index_pair > & solGIDList,
                                 const std::list< index_pair > & staGIDList) {return;}

  // Compute an estimate of the error in the integration step (BDF1-5).
  virtual double computeErrorEstimate() { return sec.ck_*ds.WRMS_errorNorm(); }

  // Interpolate solution approximation at prescribed time point (BDF1-5).
  virtual bool interpolateSolution(double timepoint,
                  N_LAS_Vector * tmpSolVectorPtr, std::vector<N_LAS_Vector*> & historyVec);

  // Interpolate MPDE solution approximation at prescribed time point (BDF1-5).
  virtual bool interpolateMPDESolution(std::vector<double>& timepoint,
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
                  const std::vector<double> & outputInterpolationTimes,
                  bool skipPrintLineOutput);

  // .SAVE output using interpolation when order is high
  virtual bool saveOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  N_LAS_Vector * solnVecPtr,
                  const double saveTime,
                  const bool doNotInterpolate); 

  // Computes the step adjustment (BDF1-5).
  virtual double computeExpoStepAdjust(double stepadjust);

  // Gets the time-integration order (BDF1-5).
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

  // 03/08/04 tscoffe:  New functions necessary for BackwardDifferentiation15
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
  //
  // Return min order of method (this should obey user option minorder)
  virtual int getMinOrder(){ return sec.minOrder_;}

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
  N_TIA_BackwardDifferentiation15(N_TIA_TIAParams & tiaP,
                                  N_TIA_StepErrorControl & secTmp,
                                  N_TIA_DataStore & dsTmp);

  // Check whether to reduce order independent of local error test
  virtual void checkReduceOrder();

  enum actionFlag { TIAAction_UNSET, TIAAction_LOWER, TIAAction_MAINTAIN, TIAAction_RAISE };

  // Used to keep track of last interpolation point in printMPDEOutputSolution
  double timept_; 
};

  // Computes the step addjustment(BDF1-5).
  // 2/16/04 tscoffe:  I'm not exactly sure what this routine is for...
inline double N_TIA_BackwardDifferentiation15::computeExpoStepAdjust(
  double stepadjust)
{
  return pow(stepadjust, 1.0 / 3.0);
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_BackwardDifferentiation15::factory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 2/16/04
//-----------------------------------------------------------------------------
inline N_TIA_TimeIntegrationMethod * N_TIA_BackwardDifferentiation15::factory
  (N_TIA_TIAParams & tiaP, 
   N_TIA_StepErrorControl & secTmp ,
   N_TIA_DataStore & dsTmp)
{
  N_TIA_TimeIntegrationMethod * timPtr =
    new N_TIA_BackwardDifferentiation15(tiaP,secTmp,dsTmp);
  return timPtr;
}

#endif     //Xyce_N_TIA_BACKWARD_DIFFERENTIATION_15_H

