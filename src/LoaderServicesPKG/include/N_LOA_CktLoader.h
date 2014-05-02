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
// Filename       : $RCSfile: N_LOA_CktLoader.h,v $
//
// Purpose        : This file contains class definitions for the loader
//                  services package.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/26/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.88 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_LOA_CktLoader_H
#define Xyce_LOA_CktLoader_H

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include <N_LOA_Loader.h>
#include <N_UTL_BreakPoint.h>
#include <N_DEV_DeviceInterface.h>


// ---------- Forward declarations --------
class N_LAS_Vector;
class N_LAS_Matrix;

//-----------------------------------------------------------------------------
// Class         : N_LOA_CktLoader
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/03/00
//-----------------------------------------------------------------------------
class N_LOA_CktLoader : public N_LOA_Loader
{
public:

  // Default constructor
  N_LOA_CktLoader();

  // Destructor
  ~N_LOA_CktLoader();

  // Method which is called to load the new-DAE contributions to
  // the nonlinear Jacobian matrix.
  bool loadDAEMatrices  (N_LAS_Vector * tmpSolVectorPtr,
                         N_LAS_Vector * tmpStaVectorPtr,
                         N_LAS_Vector * tmpStaDerivVectorPtr,
                         N_LAS_Vector * tmpStoVectorPtr,
                         N_LAS_Matrix * tmpdQdxMatrixPtr,
                         N_LAS_Matrix * tmpdFdxMatrixPtr);

  // Initializes the nonlinear problem.
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

  // Method is called to load the mask to be used in calculating error norms.
  bool loadDeviceMask();

  // Method which is called to load the new-DAE vectors, which contribute
  // to the nonlinear residual (RHS) vector.
  bool loadDAEVectors   (N_LAS_Vector * nextSolVectorPtr,
                         N_LAS_Vector * currSolVectorPtr,
                         N_LAS_Vector * lastSolVectorPtr,
                         N_LAS_Vector * nextStaVectorPtr,
                         N_LAS_Vector * currStaVectorPtr,
                         N_LAS_Vector * lastStaVectorPtr,
                         N_LAS_Vector * StaDerivVectorPtr,
                         N_LAS_Vector * nextStoVectorPtr,
                         N_LAS_Vector * currStoVectorPtr,
                         N_LAS_Vector * lastStoVectorPtr,
                         N_LAS_Vector * stoLeadCurrQCompVectorPtr,
                         N_LAS_Vector * QVectorPtr,
                         N_LAS_Vector * FVectorPtr,
                         N_LAS_Vector * dFdxdVpVectorPtr,
                         N_LAS_Vector * dQdxdVpVectorPtr);

  bool updateState      (N_LAS_Vector * nextSolVectorPtr,
                         N_LAS_Vector * currSolVectorPtr,
                         N_LAS_Vector * lastSolVectorPtr,
                         N_LAS_Vector * nextStaVectorPtr,
                         N_LAS_Vector * currStaVectorPtr,
                         N_LAS_Vector * lastStaVectorPtr,
                         N_LAS_Vector * nextStoVectorPtr,
                         N_LAS_Vector * currStoVectorPtr,
                         N_LAS_Vector * lastStoVectorPtr
                         );

  bool loadBVectorsforAC (N_LAS_Vector * bVecRealPtr,
                          N_LAS_Vector * bVecImagPtr);

  bool getBMatrixEntriesforMOR(std::vector<int>& bMatEntriesVec, std::vector<int>& bMatPosEntriesVec);

  // Function for setting the initial guess.
  bool setInitialGuess (N_LAS_Vector * solVectorPtr);

  // Function for setting a single parameter value.
    virtual bool setParam(std::string & name, double val) {
      return deviceIntPtr->setParam(name,val);
    }

    // Function for getting a single parameter value.
    virtual double getParamAndReduce(const std::string & name) {
      return deviceIntPtr->getParamAndReduce(name);
    }

    virtual bool getParamAndReduce(const std::string & name, double & val) {
      return deviceIntPtr->getParamAndReduce(name,val);
    }

  bool getVsrcLIDs (std::string & srcName, int & li_Pos, int & li_Neg, int & li_Bra);

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

  void getNumInterfaceNodes (std::vector<int> & numINodes);
  bool loadCouplingRHS(int iSubProblem, int iCouple, N_LAS_Vector * dfdvPtr);
  bool calcCouplingTerms (int iSubProblem, int iCouple, const N_LAS_Vector * dxdvPtr);
  virtual bool raiseDebugLevel (int increment);

  // Gets the time integration required breakpoint times (in a vector).
  bool getBreakPoints(std::vector< N_UTL_BreakPoint > & breakPointTimes);

  // Accessor which returns the maximum time step size (in seconds).
  double getMaxTimeStepSize();

  // Registration method for the device packaage
  bool registerDeviceInterface(N_DEV_DeviceInterface * devIntPtr);

  // Get block size from device options for block gainscale homotopy
  int getHomotopyBlockSize() const;

  // Get convergence info from devices
  bool allDevsConverged();

  // Get convergence info from inner-solves
  bool innerDevsConverged();

  // Functions needed by the NEW (power node) 2-level algorithm:
  void homotopyStepSuccess
    (const std::vector<std::string> & paramNames,
     const std::vector<double> & paramVals);

  void homotopyStepFailure ();

  void stepSuccess (int analysis);
  void stepFailure (int analysis);

  void acceptStep();

  virtual bool getInitialQnorm (std::vector<N_TIA_TwoLevelError> & tleVec );

  virtual bool getInnerLoopErrorSums (std::vector<N_TIA_TwoLevelError> & tleVec);

  bool updateStateArrays ();
  bool startTimeStep ();
  void setExternalSolverState (const N_DEV_SolverState & ss);

protected:
private :

public :

  // Pointer to the device package interface
  N_DEV_DeviceInterface * deviceIntPtr;

};

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::initializeProblem
// Purpose       : This function calls the setICs function in the device
//                 manager.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::initializeProblem(
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
  return N_LOA_CktLoader::deviceIntPtr->setICs
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

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::loadDeviceMask
// Purpose       : Load mask for use in calculating weighted norms.
// Special Notes : Loads the mask into the LAS_System.  Need only be called
//                 once.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 1/19/07
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::loadDeviceMask()
{
  return N_LOA_CktLoader::deviceIntPtr->loadDeviceMask();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::loadDAEMatrices
// Purpose       : This function
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::loadDAEMatrices  (N_LAS_Vector * tmpSolVectorPtr,
                                        N_LAS_Vector * tmpStateVectorPtr,
                                        N_LAS_Vector * tmpStateDerivVectorPtr,
                                        N_LAS_Vector * tmpStoreVectorPtr,
                                        N_LAS_Matrix * tmpdQdxMatrixPtr,
                                        N_LAS_Matrix * tmpdFdxMatrixPtr)
{
  return
    N_LOA_CktLoader::deviceIntPtr->loadDAEMatrices (tmpSolVectorPtr,
                                                    tmpStateVectorPtr,
                                                    tmpStateDerivVectorPtr,
                                                    tmpStoreVectorPtr,
                                                    tmpdQdxMatrixPtr,
                                                    tmpdFdxMatrixPtr);

}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::loadDAEVectors   (
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
                                        N_LAS_Vector * stoLeadCurrQCompVectorPtr,
                                        N_LAS_Vector * QVectorPtr,
                                        N_LAS_Vector * FVectorPtr,
                                        N_LAS_Vector * dFdxdVpVectorPtr,
                                        N_LAS_Vector * dQdxdVpVectorPtr)
{
  return N_LOA_CktLoader::deviceIntPtr->loadDAEVectors
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
       stoLeadCurrQCompVectorPtr,
       QVectorPtr,
       FVectorPtr,
       dFdxdVpVectorPtr,
       dQdxdVpVectorPtr);
}


//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 11/17/05
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::updateState
 (N_LAS_Vector * nextSolVectorPtr,
	N_LAS_Vector * currSolVectorPtr,
	N_LAS_Vector * lastSolVectorPtr,
	N_LAS_Vector * nextStaVectorPtr,
	N_LAS_Vector * currStaVectorPtr,
	N_LAS_Vector * lastStaVectorPtr,
  N_LAS_Vector * nextStoVectorPtr,
  N_LAS_Vector * currStoVectorPtr,
  N_LAS_Vector * lastStoVectorPtr
  )
{
  return deviceIntPtr->updateState
    ( nextSolVectorPtr,
      currSolVectorPtr,
      lastSolVectorPtr,
      nextStaVectorPtr,
      currStaVectorPtr,
      lastStaVectorPtr,
      nextStoVectorPtr,
      currStoVectorPtr,
      lastStoVectorPtr
      );
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::loadBVectorsforAC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::loadBVectorsforAC
  (N_LAS_Vector * bVecRealPtr,
   N_LAS_Vector * bVecImagPtr)
{
  return deviceIntPtr->loadBVectorsforAC (bVecRealPtr,
                                          bVecImagPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::getBMatrixStampforMOR
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::getBMatrixEntriesforMOR(std::vector<int>& bMatEntriesVec,
                                                     std::vector<int>& bMatPosEntriesVec)
{
  return deviceIntPtr->getBMatrixEntriesforMOR( bMatEntriesVec, bMatPosEntriesVec );
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::getInductorsEntriesforMOR
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
/*inline bool N_LOA_CktLoader::getInductorsEntriesforMOR(std::vector<int>& inductorEntriesVec)
{
  return deviceIntPtr->getInductorsEntriesforMOR( inductorEntriesVec );
}
*/
//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::setInitialGuess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::setInitialGuess(N_LAS_Vector * solVectorPtr)
{
  return N_LOA_CktLoader::deviceIntPtr->setInitialGuess (solVectorPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::getVsrcLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/05/06
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::getVsrcLIDs
  (std::string & srcName, int & li_Pos, int & li_Neg, int & li_Bra)
{
  return N_LOA_CktLoader::deviceIntPtr->getVsrcLIDs
      (srcName, li_Pos, li_Neg, li_Bra);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::updateSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::updateSources()
{
  return N_LOA_CktLoader::deviceIntPtr->updateSources();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::output
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::output ()
{
  return N_LOA_CktLoader::deviceIntPtr->output ();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::finishOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/19/04
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::finishOutput ()
{
  return N_LOA_CktLoader::deviceIntPtr->finishOutput ();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::getLinearSystemFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/17/00
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::getLinearSystemFlag()
{
  return N_LOA_CktLoader::deviceIntPtr->getLinearSystemFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::getDoubleDCOPFlag
// Purpose       : This is an accessor to allow the time integrator to determine
//                 if the current problem includes a PDE device.  If it does,
//                 then it makes sense to use two-pass DCOP calulation.  Hence
//                 the "DoubleDCOPFlag" in the name.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/25/01
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::getDoubleDCOPFlag ()
{
  return N_LOA_CktLoader::deviceIntPtr->getPDESystemFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::getLimiterFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/11/04
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::getLimiterFlag ()
{
  return N_LOA_CktLoader::deviceIntPtr->getVoltageLimiterFlag ();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/11/01
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::getBreakPoints ( std::vector<N_UTL_BreakPoint> & breakPointTimes )
{
  return N_LOA_CktLoader::deviceIntPtr->getBreakPoints(breakPointTimes);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::registerDeviceInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/08/01
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::registerDeviceInterface (N_DEV_DeviceInterface * devIntPtr)
{
  bool bsuccess = true;
  deviceIntPtr = devIntPtr;
  if (deviceIntPtr == NULL) bsuccess = false;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::getMaxTimeStepSize ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/31/01
//-----------------------------------------------------------------------------
inline double N_LOA_CktLoader::getMaxTimeStepSize ()
{
  return deviceIntPtr->getMaxTimeStepSize ();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::enablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
inline int N_LOA_CktLoader::enablePDEContinuation ()
{
  return deviceIntPtr->enablePDEContinuation();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::disablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::disablePDEContinuation ()
{
  return deviceIntPtr->disablePDEContinuation();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::getNumInterfaceNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
inline void N_LOA_CktLoader::getNumInterfaceNodes (std::vector<int> & numINodes)
{
  deviceIntPtr->getNumInterfaceNodes (numINodes);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::loadCouplingRHS
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::loadCouplingRHS(int iSubProblem, int iCouple, N_LAS_Vector * dfdvPtr)
{
  return deviceIntPtr->loadCouplingRHS (iSubProblem, iCouple, dfdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::calcCouplingTerms
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::calcCouplingTerms (int iSubProblem, int iCouple, const N_LAS_Vector * dxdvPtr)
{
  return deviceIntPtr->calcCouplingTerms (iSubProblem, iCouple, dxdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::raiseDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/23/03
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::raiseDebugLevel (int increment)
{
  return deviceIntPtr->raiseDebugLevel (increment);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::getHomotopyBlockSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL, Parallel Computational Sciences
// Creation Date : 01/26/2005
//-----------------------------------------------------------------------------
inline int N_LOA_CktLoader::getHomotopyBlockSize() const
{
  return deviceIntPtr->getHomotopyBlockSize();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::allDevsConverged
// Purpose       : Check whether any device has taken an action that renders
//                  normal convergence checks invalid (i.e. that the current
//                  step must be assumed unconverged).
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/22/05
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::allDevsConverged()
{
  return deviceIntPtr->allDevsConverged();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::innerDevsConverged
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/21/06
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::innerDevsConverged()
{
  return deviceIntPtr->innerDevsConverged();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::homotopyStepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline void N_LOA_CktLoader::homotopyStepSuccess
    (const std::vector<std::string> & paramNames,
     const std::vector<double> & paramVals)
{
  return deviceIntPtr->homotopyStepSuccess (paramNames, paramVals) ;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::homotopyStepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
inline void N_LOA_CktLoader::homotopyStepFailure ()
{
  return deviceIntPtr->homotopyStepFailure ();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
inline void N_LOA_CktLoader::stepSuccess (int analysis)
{
  return deviceIntPtr->stepSuccess (analysis);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
inline void N_LOA_CktLoader::stepFailure (int analysis)
{
  return deviceIntPtr->stepFailure (analysis);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::acceptStep
// Purpose       : Communicate to devices that the current step has been
//                 accepted
// Special Notes : Most devices need not know.  The Transmission line does.
// Scope         : public
// Creator       : Tom Russo, SNL
// Creation Date : 01/23/07
//-----------------------------------------------------------------------------
inline void N_LOA_CktLoader::acceptStep ()
{
  deviceIntPtr->acceptStep();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::getInitialQnorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/18/07
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::getInitialQnorm (std::vector<N_TIA_TwoLevelError> & tleVec )
{
  return deviceIntPtr->getInitialQnorm(tleVec);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::getInnerLoopErrorSums
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::getInnerLoopErrorSums (std::vector<N_TIA_TwoLevelError> & tleVec)
{
  return deviceIntPtr->getInnerLoopErrorSums (tleVec);
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::updateStateArrays
// Purpose       : Tells the inner loop solve to update state arrays.
// Special Notes : Needed to support voltlim with loca, with 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::updateStateArrays ()
{
  return deviceIntPtr->updateStateArrays ();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::startTimeStep ()
// Purpose       : Tells the inner loop solve to do its prediction.
// Special Notes : Needed to support 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
inline bool N_LOA_CktLoader::startTimeStep ()
{
  return deviceIntPtr->startTimeStep ();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_CktLoader::setExternalSolverState
// Purpose       :
// Special Notes : Needed to support 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/21/06
//-----------------------------------------------------------------------------
inline void N_LOA_CktLoader::setExternalSolverState (const N_DEV_SolverState & ss)
{
  return deviceIntPtr->setExternalSolverState (ss);
}

#endif
