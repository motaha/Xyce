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
// Filename      : $RCSfile: N_TIA_TimeIntegrationMethods.h,v $
//
// Purpose       : This file defines the classes for the time integration
//                 methods --- the "interface base class" along with the
//                 accompanying integration methods classes which can be
//                 used in the time integration algorithm.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.53.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_TIME_INTEG_METH_H
#define Xyce_N_TIA_TIME_INTEG_METH_H

// ---------- Standard Includes ----------
#include <iostream>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <list>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_UTL_Xyce.h>
#include <N_IO_fwd.h>
#include <N_TIA_DataStore.h>
#include <N_ANP_SweepParam.h>
#include <N_TIA_TimeIntegrationEnum.h>

#include <N_ANP_OutputMgrAdapter.h>

// ----------   Forward Declarations ---------
class N_PDS_Comm;
class N_TIA_StepErrorControl;
class N_TIA_TimeIntInfo;
class N_LAS_Vector;
class N_LAS_Matrix;
class N_ERH_ErrorMgr;

//-----------------------------------------------------------------------------
// Class         : N_TIA_TimeIntegrationMethod
//
// Purpose       : This is the integration methods base class, from which
//                 specific integration methods (such as BDF15, trap, etc) are
//                 derrived.
//
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------

class N_TIA_TimeIntegrationMethod
{

public:

  // Default constructor.
  N_TIA_TimeIntegrationMethod( N_TIA_TIAParams & tiaP,
                               N_TIA_StepErrorControl & secTmp,
                               N_TIA_DataStore & dsTmp);

  // Destructor
  virtual ~N_TIA_TimeIntegrationMethod();

  // Predict solution at next time point (abstract).
  virtual void obtainPredictor() = 0;

  // Evaluate the predictor derivative formula (abstract).
  virtual void obtainPredictorDeriv() = 0;

  // Evaluate the corrector derivative formula (abstract).
  virtual void obtainCorrectorDeriv() = 0;

  // Compute an estimate of the error in the integration step (abstract).
  virtual void updateDerivsBlock(const list < index_pair > & solGIDList,
                                 const list < index_pair > & staGIDList) = 0;

  // Compute an estimate of the error in the integration step (abstract).
  virtual double computeErrorEstimate() = 0;

  // Interpolate solution, state or store approximation at prescribed time point (abstract).
  virtual bool interpolateSolution(double timepoint,
      	N_LAS_Vector * tmpSolVectorPtr, vector<N_LAS_Vector*> & historyVec) = 0;

  // Print output using interpolation when order is high
  virtual bool printOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  const double time,
                  N_LAS_Vector * solnVecPtr,
                  const bool doNotInterpolate,
                  const vector<double> &outputInterpolationTimes,
                  bool skipPrintLineOutput );

  // Print MPDE output using local interpolation methods
  virtual bool printMPDEOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  const double time,
                  N_LAS_Vector * solnVecPtr,
                  const std::vector<double>& fastTimes );

  // Print WaMPDE output using local interpolation methods
  virtual bool printWaMPDEOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  const double time,
                  N_LAS_Vector * solnVecPtr,
                  const std::vector<double>& fastTimes,
                  const int phiGID );

  // .SAVE output using interpolation when order is high
  virtual bool saveOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  N_LAS_Vector * solnVecPtr,
                  const double saveTime,
                  const bool doNotInterpolate);

  // Computes the step addjustment.
  virtual double computeExpoStepAdjust(double stepadjust) = 0;

  // Gets the time-integration order (abstract).
  virtual int getOrder() = 0;
  virtual int getNumberOfSteps () { return 0; }
  virtual int getUsedOrder () = 0;
  virtual int getNscsco () { return 0; }

  virtual void getInitialQnorm (N_TIA_TwoLevelError & tle) = 0;
  virtual void setupTwoLevelError(N_TIA_TwoLevelError & tle) = 0;

  // This is for new-DAE.
  virtual void updateStateDeriv () {}

  // calculates dQ/dt component of store vector and adds it to store vector
  virtual void updateLeadCurrent () {}

  // Returns the current partial time derivative for either the solution or
  // state vector.
  virtual double partialTimeDeriv();

  // Gets the leading coefficient for the specified time-integration method.
  virtual double getLeadingCoeff() { return leadingCoeff; }

  // sets the leading coefficient for the specified time-integration method.
  virtual void setLeadingCoeff(double & LC) { leadingCoeff = LC; }

  // 03/08/04 tscoffe:  New functions necessary for BackwardDifferentiation15
  // Evaluate residual for nonlinear solver
  virtual void obtainResidual();

  // Evaluate Jacobian for nonlinear solver
  virtual void obtainJacobian();

  // Apply Jacobian for nonlinear solver
  virtual void applyJacobian(const N_LAS_Vector& input, N_LAS_Vector& result);

  // Update history array after a successful step
  virtual void updateHistory();

  // Restore history array after a failed step
  virtual void restoreHistory();

  // Return max order of method (this should obey user option maxorder)
  virtual int getMaxOrder();

  // Update method coefficients
  virtual void updateCoeffs();

  // Initialize method with initial solution & step-size
  virtual void initialize();

  // setup 2-level data.
  virtual void setTwoLevelTimeInfo(const N_TIA_TimeIntInfo & tiInfo);

  // Reject a step (this turns back history and chooses new order & step-size)
  virtual void rejectStep();

  virtual void rejectStepForHabanero () {};

  // Complete a step (this updates history and chooses new order & step-size)
  virtual void completeStep();

protected:
private :

public :

  // Reference to the TIA data-store object.
  N_TIA_DataStore & ds;

  // Reference to step-error control object.
  N_TIA_StepErrorControl & sec;

  // Time-integration method leading coefficient value.
  double leadingCoeff;

  // Reference to TIA params object.
  N_TIA_TIAParams & tiaParams;
};

//-----------------------------------------------------------------------------
// Class         : N_TIA_WorkingIntegrationMethod
// Purpose       : This class provides a way for obtaining a specific
//                 working integration method and its associated data items.
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class N_TIA_WorkingIntegrationMethod
{
  public:

    // Default constructor.
    N_TIA_WorkingIntegrationMethod(N_TIA_TIAParams & tiaP,
                                  N_TIA_StepErrorControl & secTmp,
                                  N_TIA_DataStore & dsTmp);

    // Constructor
    N_TIA_WorkingIntegrationMethod(
    const unsigned int integMethod, N_TIA_TIAParams & tiaP,
                                    N_TIA_StepErrorControl & secTmp,
                                    N_TIA_DataStore & dsTmp);

    // Destructor
    ~N_TIA_WorkingIntegrationMethod();

    // Method which creates a time-integration method.
    void createTimeIntegMethod(const unsigned int integMethod);

#ifdef Xyce_VERBOSE_TIME
    // Output the current time-integration method.
    void printWorkingIntegMethod();
#endif

    // Current integration method flag.
    unsigned int workingIntegMethod;

    // accessors to the integration method object functions:
    N_TIA_TimeIntegrationMethod * getIntegMethodPtr() { return integMethodPtr; }

    double partialTimeDeriv() { return integMethodPtr->partialTimeDeriv(); }
    void obtainPredictor() { integMethodPtr->obtainPredictor(); }
    void obtainPredictorDeriv() { integMethodPtr->obtainPredictorDeriv(); }
    void obtainCorrectorDeriv() { integMethodPtr->obtainCorrectorDeriv(); }
    void updateDerivsBlock ( const list<index_pair> & solGIDList, const list<index_pair> & staGIDList)
    { integMethodPtr->updateDerivsBlock (solGIDList, staGIDList); }

    int getOrder() { return integMethodPtr->getOrder(); }
    int getUsedOrder() { return integMethodPtr->getUsedOrder(); }
    int getNumberOfSteps() { return integMethodPtr->getNumberOfSteps(); }
    int getNscsco () { return integMethodPtr->getNscsco(); }
    void getInitialQnorm (N_TIA_TwoLevelError & tle) { return integMethodPtr->getInitialQnorm (tle); }
    void setupTwoLevelError(N_TIA_TwoLevelError & tle) { return integMethodPtr->setupTwoLevelError(tle); }
    void setTwoLevelTimeInfo(const N_TIA_TimeIntInfo & tiInfo) { return integMethodPtr->setTwoLevelTimeInfo(tiInfo); }
    void updateCoeffs() { return integMethodPtr->updateCoeffs(); }
    void rejectStepForHabanero () { return integMethodPtr->rejectStepForHabanero(); }
    void initialize() { return integMethodPtr->initialize(); }
    void completeStep() { return integMethodPtr->completeStep(); }
    void rejectStep() { return integMethodPtr->rejectStep(); }
    double computeErrorEstimate() { return integMethodPtr->computeErrorEstimate(); }
    void updateStateDeriv () { return integMethodPtr->updateStateDeriv (); }
    void updateLeadCurrent () { return integMethodPtr->updateLeadCurrent (); }
    void obtainResidual() { return integMethodPtr->obtainResidual(); }
    void obtainJacobian() { return integMethodPtr->obtainJacobian(); }
    void applyJacobian(const N_LAS_Vector& input, N_LAS_Vector& result)
    { return integMethodPtr->applyJacobian(input, result); }

    bool printMPDEOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  const double time,
                  N_LAS_Vector * solnVecPtr,
                  const std::vector<double>& fastTimes )
    {
      return integMethodPtr->printMPDEOutputSolution(
                  outputMgrAdapterRCPtr, time, solnVecPtr, fastTimes );
    }

    bool printWaMPDEOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  const double time,
                  N_LAS_Vector * solnVecPtr,
                  const std::vector<double>& fastTimes,
                  const int phiGID )
    {
      return integMethodPtr->printWaMPDEOutputSolution(
                  outputMgrAdapterRCPtr, time, solnVecPtr, fastTimes, phiGID );
    }

    bool printOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  const double time,
                  N_LAS_Vector * solnVecPtr,
                  const bool doNotInterpolate,
                  const vector<double> &outputInterpolationTimes,
                  bool skipPrintLineOutput )
    {
      return integMethodPtr->printOutputSolution(
          outputMgrAdapterRCPtr, time, solnVecPtr, doNotInterpolate, outputInterpolationTimes, skipPrintLineOutput) ;
    }

    bool saveOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  N_LAS_Vector * solnVecPtr,
                  const double saveTime,
                  const bool doNotInterpolate)
    {
      return integMethodPtr->saveOutputSolution(
          outputMgrAdapterRCPtr, solnVecPtr, saveTime, doNotInterpolate);
    }

    // TIA params object.
    N_TIA_TIAParams & tiaParams;

    // Reference to the TIA data-store object.
    N_TIA_DataStore & ds;

    // Reference to step-error control object.
    N_TIA_StepErrorControl & sec;

  private:
    // Pointer to the integration method.
    N_TIA_TimeIntegrationMethod * integMethodPtr;
};




//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::obtainJacobian
// Purpose       : Evaluate Jacobian for nonlinear solver
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 3/08/04
//-----------------------------------------------------------------------------
inline void N_TIA_TimeIntegrationMethod::obtainJacobian()
{
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::applyJacobian
// Purpose       : Apply Jacobian for nonlinear solver
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/29/08
//-----------------------------------------------------------------------------
inline void N_TIA_TimeIntegrationMethod::applyJacobian(const N_LAS_Vector& input, N_LAS_Vector& result)
{
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::updateHistory
// Purpose       : Update history array after a successful step
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 3/08/04
//-----------------------------------------------------------------------------
inline void N_TIA_TimeIntegrationMethod::updateHistory()
{
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::restoreHistory
// Purpose       : Restore history array after a failed step
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 3/08/04
//-----------------------------------------------------------------------------
inline void N_TIA_TimeIntegrationMethod::restoreHistory()
{
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::getMaxOrder
// Purpose       : Return max order of method (this should obey user option maxorder)
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 3/08/04
//-----------------------------------------------------------------------------
inline int N_TIA_TimeIntegrationMethod::getMaxOrder()
{
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::updateCoeffs
// Purpose       : Update method coefficients
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 3/08/04
//-----------------------------------------------------------------------------
inline void N_TIA_TimeIntegrationMethod::updateCoeffs()
{
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::initialize
// Purpose       : Initialize method with initial solution & step-size
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 3/09/04
//-----------------------------------------------------------------------------
inline void N_TIA_TimeIntegrationMethod::initialize()
{
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::setTwoLevelTimeInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/01/07
//-----------------------------------------------------------------------------
inline void N_TIA_TimeIntegrationMethod::setTwoLevelTimeInfo(const N_TIA_TimeIntInfo & tiInfo)
{
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::printMPDEOutputSolution()
// Purpose       : Use local time integration algorithm interpolants to
//                 print transient MPDE output.
// Scope         : public
// Creator       : Todd Coffey, SNL, 1414
// Creation Date : 12/05/06
//-----------------------------------------------------------------------------
inline bool N_TIA_TimeIntegrationMethod::printMPDEOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  const double time,
                  N_LAS_Vector * solnVecPtr,
                  const std::vector<double>& fastTimes )
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_TimeIntegrationMethod::printWaMPDEOutputSolution()
// Purpose       : Use local time integration algorithm interpolants to
//                 print transient WaMPDE output.
// Scope         : public
// Creator       : Todd Coffey, SNL, 1414
// Creation Date : 12/15/06
//-----------------------------------------------------------------------------
inline bool N_TIA_TimeIntegrationMethod::printWaMPDEOutputSolution(
                  RefCountPtr< N_ANP_OutputMgrAdapter > outputMgrAdapterRCPtr,
                  const double time,
                  N_LAS_Vector * solnVecPtr,
                  const std::vector<double>& fastTimes,
                  const int phiGID )
{
  return false;
}

#endif // Xyce_N_TIA_TIME_INTEG_METH_H

