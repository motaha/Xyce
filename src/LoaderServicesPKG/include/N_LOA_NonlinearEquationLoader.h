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
// Filename       : $RCSfile: N_LOA_NonlinearEquationLoader.h,v $
//
// Purpose        : This file contains the interface class that sits between
//                  the nonlinear solver and (ultimately) the time integrator.
//                  For transient calculations only the time integrator knows
//                  what to sum into the Jacobian and Residual vector.
//
// Special Notes  :
//
// Creator        : Todd Coffey, SNL
//
// Creation Date  : 07/29/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.13.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:46 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_LOA_NonlinearEquationLoader_H
#define Xyce_LOA_NonlinearEquationLoader_H

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>
#include <N_LOA_Loader.h>
#include <N_UTL_BreakPoint.h>

#include <N_DEV_DeviceInterface.h>
#include <N_ANP_AnalysisInterface.h>

// ---------- Forward declarations --------
class N_ANP_AnalysisInterface;
class N_LAS_Vector;
class N_LAS_Matrix;

//-----------------------------------------------------------------------------
// Class         : N_LOA_NonlinearEquationLoader
//
// Purpose        : This class contains the interface class that sits between
//                  the nonlinear solver and (ultimately) the time integrator.
//                  For transient calculations only the time integrator knows
//                  what to sum into the Jacobian and Residual vector.
//
// Special Notes :
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
class N_LOA_NonlinearEquationLoader : public N_LOA_Loader
{
public:

  // Default constructor
  N_LOA_NonlinearEquationLoader();

  // Destructor
  ~N_LOA_NonlinearEquationLoader();

  // Method which is called to load the nonlinear Jacobian matrix.
  bool loadJacobian ();

  // Method which is called to apply the nonlinear Jacobian matrix.
  bool applyJacobian (const N_LAS_Vector& input, N_LAS_Vector& result);

#if 0
  // Initializes the nonlinear problem.
  bool initializeProblem();
#else
  bool initializeProblem( N_LAS_Vector * nextSolVectorPtr,
                          N_LAS_Vector * currSolVectorPtr,
                          N_LAS_Vector * lastSolVectorPtr,
                          N_LAS_Vector * nextStaVectorPtr,
                          N_LAS_Vector * currStaVectorPtr,
                          N_LAS_Vector * lastStaVectorPtr,
                          N_LAS_Vector * StateDerivVectorPtr,
                          N_LAS_Vector * nextStoVectorPtr,
                          N_LAS_Vector * currStoVectorPtr,
                          N_LAS_Vector * lastStoVectorPtr,
                          N_LAS_Vector * QVectorPtr,
                          N_LAS_Vector * FVectorPtr,
                          N_LAS_Vector * dFdxdVpVectorPtr,
                          N_LAS_Vector * dQdxdVpVectorPtr);
#endif

  // Method which is called to load the nonlinear residual (RHS) vector.
  bool loadRHS ();

  // Function for setting the initial guess.
  bool setInitialGuess (N_LAS_Vector * solVectorPtr);

  // Function for setting a single parameter value.
  bool setParam (string & name, double val);

  // Function for getting a single parameter value.
  double getParam (const string & name);
  bool getParam (const string & name, double & val);
  bool getVsrcLIDs (string & srcName, int & li_Pos, int & li_Neg, int & li_Bra);

  // Method which is called to update the sources.
  bool updateSources();
  bool getLinearSystemFlag();

  // Get the voltage limiter flag:
  bool getLimiterFlag ();

  // Gets the double DC Operating Point flag - used for PDE devices.
  bool getDoubleDCOPFlag();
  bool output();
  bool finishOutput();

  // two-level newton functions:
  int  enablePDEContinuation ();
  bool disablePDEContinuation ();

  void getNumInterfaceNodes (vector<int> & numINodes);
  bool loadCouplingRHS   (int iSubProblem, int iCouple, N_LAS_Vector * dfdvPtr);
  bool calcCouplingTerms (int iSubProblem, int iCouple, const N_LAS_Vector * dxdvPtr);
  virtual bool raiseDebugLevel (int increment);

  // Gets the time integration required breakpoint times (in a vector).
  bool getBreakPoints(vector < N_UTL_BreakPoint > & breakPointTimes);

  // Accessor which returns the maximum time step size (in seconds).
  double getMaxTimeStepSize();

  // Accessors for timing information

  // Gets the nonlinear residual load time.
  double getResidualTime();

  // Gets the nonlinear Jacobian load time.
  double getJacobianTime();

  // Registration method for the device packaage
  bool registerDeviceInterface(N_DEV_DeviceInterface * devIntPtr);

  // Registration method for the device packaage
  bool registerTIA (N_ANP_AnalysisInterface * tiaPtr);

  // Get block size from device options for block gainscale homotopy
  int getHomotopyBlockSize() const;

  // Get convergence info from devices
  bool allDevsConverged();

  // Get convergence info from inner-solves
  bool innerDevsConverged();

protected:
private :

public :

  // Pointer to the device package interface
  N_DEV_DeviceInterface * deviceIntPtr;

  // Pointer to the time integration manager.
  N_ANP_AnalysisInterface * tiaPtr;
};

#if 0
//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::initializeProblem
// Purpose       : This function calls the setICs function in the device
//                 manager.
//
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::initializeProblem()
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->setICs();
}
#else
inline bool N_LOA_NonlinearEquationLoader::initializeProblem(
                          N_LAS_Vector * nextSolVectorPtr,
                          N_LAS_Vector * currSolVectorPtr,
                          N_LAS_Vector * lastSolVectorPtr,
                          N_LAS_Vector * nextStaVectorPtr,
                          N_LAS_Vector * currStaVectorPtr,
                          N_LAS_Vector * lastStaVectorPtr,
                          N_LAS_Vector * StateDerivVectorPtr,
                          N_LAS_Vector * nextStoVectorPtr,
                          N_LAS_Vector * currStoVectorPtr,
                          N_LAS_Vector * lastStoVectorPtr,
                          N_LAS_Vector * QVectorPtr,
                          N_LAS_Vector * FVectorPtr,
                          N_LAS_Vector * dFdxdVpVectorPtr,
                          N_LAS_Vector * dQdxdVpVectorPtr)
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->setICs
      (nextSolVectorPtr,
       currSolVectorPtr,
       lastSolVectorPtr,
       nextStaVectorPtr,
       currStaVectorPtr,
       lastStaVectorPtr,
       StateDerivVectorPtr,
       nextStoVectorPtr,
       currStoVectorPtr,
       lastStoVectorPtr,
       QVectorPtr,
       FVectorPtr,
       dFdxdVpVectorPtr,
       dQdxdVpVectorPtr);
}
#endif

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::loadJacobian
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::loadJacobian ()
{
  return (N_LOA_NonlinearEquationLoader::tiaPtr->loadJacobian ());
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::applyJacobian
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::applyJacobian (const N_LAS_Vector& input, N_LAS_Vector& result)
{
  return (N_LOA_NonlinearEquationLoader::tiaPtr->applyJacobian(input,result));
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::loadRHS
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::loadRHS ()
{
  return (N_LOA_NonlinearEquationLoader::tiaPtr->loadRHS());
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::setInitialGuess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::setInitialGuess(N_LAS_Vector * solVectorPtr)
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->setInitialGuess (solVectorPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::setParam
// Purpose       :
// Special Notes : Used for continuation calculations.
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::setParam(string & name, double val)
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->setParam(name,val);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::getParam
// Purpose       :
// Special Notes : Used for continuation calculations.
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline double N_LOA_NonlinearEquationLoader::getParam (const string & name)
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->getParam(name);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::getParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::getParam (const string & name, double & val)
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->getParam(name,val);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::getVsrcLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::getVsrcLIDs
  (string & srcName, int & li_Pos, int & li_Neg, int & li_Bra)
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->getVsrcLIDs
        (srcName, li_Pos, li_Neg, li_Bra);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::updateSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::updateSources()
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->updateSources();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::output
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::output ()
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->output ();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::finishOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::finishOutput ()
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->finishOutput ();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::getLinearSystemFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::getLinearSystemFlag()
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->getLinearSystemFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::getDoubleDCOPFlag
// Purpose       : This is an accessor to allow the time integrator to determine
//                 if the current problem includes a PDE device.  If it does,
//                 then it makes sense to use two-pass DCOP calulation.  Hence
//                 the "DoubleDCOPFlag" in the name.
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::getDoubleDCOPFlag ()
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->getPDESystemFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::getLimiterFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::getLimiterFlag ()
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->getVoltageLimiterFlag ();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::getResidualTime
// Purpose       : Returns timing information from residual load.
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline double N_LOA_NonlinearEquationLoader::getResidualTime()
{
  return N_LOA_NonlinearEquationLoader::tiaPtr->getResidualTime();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::getJacobianTime
// Purpose       : Returns timing information from Jacobian load.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline double N_LOA_NonlinearEquationLoader::getJacobianTime()
{
  return N_LOA_NonlinearEquationLoader::tiaPtr->getJacobianTime();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::getBreakPoints ( vector<N_UTL_BreakPoint> & breakPointTimes )
{
  return N_LOA_NonlinearEquationLoader::deviceIntPtr->getBreakPoints(breakPointTimes);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::registerDeviceInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::registerDeviceInterface (N_DEV_DeviceInterface * devIntPtr)
{
  bool bsuccess = true;

  deviceIntPtr = devIntPtr;

  if (deviceIntPtr == NULL) bsuccess = false;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::getMaxTimeStepSize ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline double N_LOA_NonlinearEquationLoader::getMaxTimeStepSize ()
{
  return deviceIntPtr->getMaxTimeStepSize ();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::enablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline int N_LOA_NonlinearEquationLoader::enablePDEContinuation ()
{
  return deviceIntPtr->enablePDEContinuation();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::disablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::disablePDEContinuation ()
{
  return deviceIntPtr->disablePDEContinuation();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::getNumInterfaceNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline void N_LOA_NonlinearEquationLoader::getNumInterfaceNodes (vector<int> & numINodes)
{
  deviceIntPtr->getNumInterfaceNodes (numINodes);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::loadCouplingRHS
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::loadCouplingRHS(int iSubProblem, int iCouple, N_LAS_Vector * dfdvPtr)
{
  return deviceIntPtr->loadCouplingRHS (iSubProblem, iCouple, dfdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::calcCouplingTerms
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::calcCouplingTerms (int iSubProblem, int iCouple, const N_LAS_Vector * dxdvPtr)
{
  return deviceIntPtr->calcCouplingTerms (iSubProblem, iCouple, dxdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::raiseDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::raiseDebugLevel (int increment)
{
  return deviceIntPtr->raiseDebugLevel (increment);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::getHomotopyBlockSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline int N_LOA_NonlinearEquationLoader::getHomotopyBlockSize() const
{
  return deviceIntPtr->getHomotopyBlockSize();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::allDevsConverged
// Purpose       : Check whether any device has taken an action that renders
//                  normal convergence checks invalid (i.e. that the current
//                  step must be assumed unconverged).
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::allDevsConverged()
{
  return deviceIntPtr->allDevsConverged();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::innerDevsConverged
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::innerDevsConverged()
{
  return deviceIntPtr->innerDevsConverged();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::registerTIA
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
inline bool N_LOA_NonlinearEquationLoader::registerTIA  (N_ANP_AnalysisInterface * tmp_tiaPtr)
{
  bool bsuccess = true;

  N_LOA_NonlinearEquationLoader::tiaPtr = tmp_tiaPtr;

  if (N_LOA_NonlinearEquationLoader::tiaPtr == NULL) bsuccess = false;

  return bsuccess;
}

#endif

