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
// Filename      : $RCSfile: N_TIA_DataStore.C,v $
//
// Purpose       : This file creates & initializes the data arrays needed for
//                 the time integration algorithms.
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
// Revision Number: $Revision: 1.136.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

#include <iostream>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

// ----------   Xyce Includes   ----------

#include <N_TIA_DataStore.h>
#include <N_TIA_TIAParams.h>

#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>

#include <N_ERH_ErrorMgr.h>
#include <Teuchos_Utils.hpp>

// ---------- Forward Declarations ----------

// ---------- Static Initializations ----------

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::N_TIA_DataStore
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------

N_TIA_DataStore::N_TIA_DataStore(N_TIA_TIAParams * tiaPtr, N_LAS_System * lsPtr)
  :
  limiterFlag(false),
  lasSysPtr(lsPtr),
  solutionSize(tiaPtr->solutionSize),
  stateSize(tiaPtr->stateSize),
  tiaParamsPtr_(tiaPtr),
  nextSolPtrSwitched(false),
  tmpSolVectorPtr(0),
  tmpStaVectorPtr(0),
  tmpStaDerivPtr(0),
  tmpStaDivDiffPtr(0),
  tmpStoVectorPtr(0),
  xn0Ptr (0),
  currSolutionPtr(0),
  lastSolutionPtr(0),
  oldeSolutionPtr(0),
  nextSolutionPtr(0),
  flagSolutionPtr(0),
  savedNextSolutionPtr(0),
  currStatePtr(0),
  lastStatePtr(0),
  oldeStatePtr(0),
  nextStatePtr(0),
  currStorePtr(0),
  lastStorePtr(0),
  oldeStorePtr(0),
  nextStorePtr(0),
  currStoreLeadCurrQCompPtr(0),
  lastStoreLeadCurrQCompPtr(0),
  oldStoreLeadCurrQCompPtr(0),
  nextStoreLeadCurrQCompPtr(0),
  currStoreLeadCurrQCompDerivPtr(0),
  lastStoreLeadCurrQCompDerivPtr(0),
  oldStoreLeadCurrQCompDerivPtr(0),
  nextStoreLeadCurrQCompDerivPtr(0),
  currSolutionDerivPtr(0),
  lastSolutionDerivPtr(0),
  oldeSolutionDerivPtr(0),
  nextSolutionDerivPtr(0),
  currStateDerivPtr(0),
  lastStateDerivPtr(0),
  oldeStateDerivPtr(0),
  nextStateDerivPtr(0),
  currSolutionDivDiffPtr(0),
  lastSolutionDivDiffPtr(0),
  oldeSolutionDivDiffPtr(0),
  nextSolutionDivDiffPtr(0),
  currStateDivDiffPtr(0),
  lastStateDivDiffPtr(0),
  oldeStateDivDiffPtr(0),
  nextStateDivDiffPtr(0),
  errWtVecPtr(0),
  absErrTolPtr(0),
  relErrTolPtr(0),
  JMatrixPtr(0),
  RHSVectorPtr(0),
#ifdef Xyce_DEBUG_DEVICE
  JdxpVectorPtr(0),
#endif
  newtonCorrectionPtr(0),
  deviceMaskPtr(0),
  indexVecsInitialized(false),
 
  qErrWtVecPtr(0),
  daeQVectorPtr(0),
  daeFVectorPtr(0),
  dQdxMatrixPtr(0),
  dFdxMatrixPtr(0),
  dQdxVecVectorPtr(0),
  dFdxVecVectorPtr(0),
  dFdxdVpVectorPtr(0),
  dQdxdVpVectorPtr(0),
  qn0Ptr(0),
  qpn0Ptr(0),
  sn0Ptr(0),
  spn0Ptr(0),
  ston0Ptr(0),
  stopn0Ptr(0),
  qNewtonCorrectionPtr(0),
  sNewtonCorrectionPtr(0),
  stoNewtonCorrectionPtr(0),
  stoLeadCurrQCompNewtonCorrectionPtr(0),
  delta_x(0),
  delta_q(0),
  tmpXn0APtr(0),
  tmpXn0BPtr(0)
{
  // temporary vectors:
  tmpSolVectorPtr = lasSysPtr->builder().createVector();
  tmpStaVectorPtr = lasSysPtr->builder().createStateVector();
  tmpStaDerivPtr = lasSysPtr->builder().createStateVector();
  tmpStaDivDiffPtr = lasSysPtr->builder().createStateVector();
  tmpStoVectorPtr = lasSysPtr->builder().createStoreVector();

  xn0Ptr = lasSysPtr->builder().createVector();

  // solution vectors:
  currSolutionPtr = lasSysPtr->builder().createVector();
  lastSolutionPtr = lasSysPtr->builder().createVector();
  nextSolutionPtr = lasSysPtr->builder().createVector();
  flagSolutionPtr = lasSysPtr->builder().createVector();

  // state vectors:
  currStatePtr    = lasSysPtr->builder().createStateVector();
  lastStatePtr    = lasSysPtr->builder().createStateVector();
  nextStatePtr    = lasSysPtr->builder().createStateVector();

  // store vectors:
  currStorePtr    = lasSysPtr->builder().createStoreVector();
  lastStorePtr    = lasSysPtr->builder().createStoreVector();
  nextStorePtr    = lasSysPtr->builder().createStoreVector();
  currStoreLeadCurrQCompPtr = lasSysPtr->builder().createStoreVector();
  lastStoreLeadCurrQCompPtr = lasSysPtr->builder().createStoreVector();
  //oldStoreLeadCurrQCompPtr = lasSysPtr->builder().createStoreVector();
  nextStoreLeadCurrQCompPtr = lasSysPtr->builder().createStoreVector();

  // solution derivative vectors:
  currSolutionDerivPtr = lasSysPtr->builder().createVector();
  lastSolutionDerivPtr = lasSysPtr->builder().createVector();
  nextSolutionDerivPtr = lasSysPtr->builder().createVector();

  // state derivative vectors:
  currStateDerivPtr = lasSysPtr->builder().createStateVector();
  lastStateDerivPtr = lasSysPtr->builder().createStateVector();
  nextStateDerivPtr = lasSysPtr->builder().createStateVector();
  
  // store derivative vec for lead currents.
  currStoreLeadCurrQCompDerivPtr = lasSysPtr->builder().createStoreVector();
  lastStoreLeadCurrQCompDerivPtr = lasSysPtr->builder().createStoreVector();
  //oldStoreLeadCurrQCompDerivPtr = lasSysPtr->builder().createStoreVector();
  nextStoreLeadCurrQCompDerivPtr = lasSysPtr->builder().createStoreVector();

  // solution scaled divided differences:
  currSolutionDivDiffPtr = lasSysPtr->builder().createVector();
  lastSolutionDivDiffPtr = lasSysPtr->builder().createVector();
  nextSolutionDivDiffPtr = lasSysPtr->builder().createVector();

  // state scaled divided differences:
  currStateDivDiffPtr = lasSysPtr->builder().createStateVector();
  lastStateDivDiffPtr = lasSysPtr->builder().createStateVector();
  nextStateDivDiffPtr = lasSysPtr->builder().createStateVector();

  // error vectors:
  errWtVecPtr      = lasSysPtr->builder().createVector();
  absErrTolPtr     = lasSysPtr->builder().createVector();
  relErrTolPtr     = lasSysPtr->builder().createVector();

  errWtVecPtr->putScalar(1.0);

  deviceMaskPtr     = lasSysPtr->builder().createVector();
  deviceMaskPtr->putScalar(1.0);

  absErrTolPtr->putScalar(tiaParamsPtr_->absErrorTol);
  relErrTolPtr->putScalar(tiaParamsPtr_->relErrorTol);

  // nonlinear solution vectors:
  newtonCorrectionPtr = lasSysPtr->builder().createVector();
  
  // new-DAE stuff:
  // Error Vectors
  qErrWtVecPtr      = lasSysPtr->builder().createVector();

  // DAE formulation vectors
  daeQVectorPtr      = lasSysPtr->builder().createVector();
  daeFVectorPtr      = lasSysPtr->builder().createVector();

  // DAE formulation matrices
  dQdxMatrixPtr = lasSysPtr->builder().createMatrix();
  dFdxMatrixPtr = lasSysPtr->builder().createMatrix();

  dQdxVecVectorPtr = lasSysPtr->builder().createVector();
  dFdxVecVectorPtr = lasSysPtr->builder().createVector();

  // History arrays 

//  if tiaParamsPtr_->method == 7 
  int sizeOfHistory = tiaParamsPtr_->maxOrder+1;
//  else
//    int sizeOfHistory = 3;

  for (int i=0;i<sizeOfHistory;++i)
  {
    xHistory.push_back(lasSysPtr->builder().createVector());
    qHistory.push_back(lasSysPtr->builder().createVector());
    sHistory.push_back(lasSysPtr->builder().createStateVector());
    stoHistory.push_back(lasSysPtr->builder().createStoreVector());
    stoLeadCurrQCompHistory.push_back(lasSysPtr->builder().createStoreVector());
  }
//  N_LAS_Vector * tmpVectorPtr = lasSysPtr->builder().createVector();
//  xHistory.resize(sizeOfHistory,*tmpVectorPtr);
//  qHistory.resize(sizeOfHistory,*tmpVectorPtr);
//  N_LAS_Vector * tmpStaVectorPtr = lasSysPtr->builder().createStateVector();
//  sHistory.resize(sizeOfHistory,*tmpStaVectorPtr);
//  delete tmpVectorPtr;
//  delete tmpStaVectorPtr;

  // Predictors
  qn0Ptr = lasSysPtr->builder().createVector();
  qpn0Ptr = lasSysPtr->builder().createVector();
  sn0Ptr = lasSysPtr->builder().createStateVector();
  spn0Ptr = lasSysPtr->builder().createStateVector();
  ston0Ptr = lasSysPtr->builder().createStoreVector();
  stopn0Ptr = lasSysPtr->builder().createStoreVector();
  stoQCn0Ptr = lasSysPtr->builder().createStoreVector();
  stoQCpn0Ptr = lasSysPtr->builder().createStoreVector();

  // Nonlinear solution vector:
  qNewtonCorrectionPtr = lasSysPtr->builder().createVector();
  sNewtonCorrectionPtr = lasSysPtr->builder().createStateVector();
  stoNewtonCorrectionPtr = lasSysPtr->builder().createStoreVector();
  stoLeadCurrQCompNewtonCorrectionPtr = lasSysPtr->builder().createStoreVector();

  // Step-size selection temporary vectors
  delta_x = lasSysPtr->builder().createVector();
  delta_q = lasSysPtr->builder().createVector();

  // Temporary vector for MPDE & WaMPDE interpolation
  tmpXn0APtr = lasSysPtr->builder().createVector();
  tmpXn0BPtr = lasSysPtr->builder().createVector();

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::N_TIA_DataStore
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------

N_TIA_DataStore::~N_TIA_DataStore()

{
  if (tmpSolVectorPtr) { delete tmpSolVectorPtr; tmpSolVectorPtr=0;}
  if (tmpStaVectorPtr) { delete tmpStaVectorPtr; tmpStaVectorPtr=0;}
  if (tmpStaDerivPtr) { delete tmpStaDerivPtr;   tmpStaDerivPtr=0;}
  if (tmpStaDivDiffPtr) { delete tmpStaDivDiffPtr; tmpStaDivDiffPtr=0;}
  if (tmpStoVectorPtr) { delete tmpStoVectorPtr; tmpStoVectorPtr=0;}
  if (xn0Ptr) { delete xn0Ptr; xn0Ptr=0;}

  if (currSolutionPtr) { delete currSolutionPtr; currSolutionPtr=0;}
  if (lastSolutionPtr) { delete lastSolutionPtr; lastSolutionPtr=0;}

  if (nextSolutionPtr) { delete nextSolutionPtr; nextSolutionPtr=0;}
  if (flagSolutionPtr) { delete flagSolutionPtr; flagSolutionPtr=0;}

  if (currStatePtr) { delete currStatePtr; currStatePtr=0;}
  if (lastStatePtr) { delete lastStatePtr; lastStatePtr=0;}

  if (nextStatePtr) { delete nextStatePtr; nextStatePtr=0;}

  if (currStorePtr) { delete currStorePtr; currStorePtr=0;}
  if (lastStorePtr) { delete lastStorePtr; lastStorePtr=0;}
  if (nextStorePtr) { delete nextStorePtr; nextStorePtr=0;}
  if (currStoreLeadCurrQCompPtr) { delete currStoreLeadCurrQCompPtr; currStoreLeadCurrQCompPtr=0;}
  if (lastStoreLeadCurrQCompPtr) { delete lastStoreLeadCurrQCompPtr; lastStoreLeadCurrQCompPtr=0;}
  //if (oldStoreLeadCurrQCompPtr) { delete oldStoreLeadCurrQCompPtr; oldStoreLeadCurrQCompPtr=0;}
  if (nextStoreLeadCurrQCompPtr) { delete nextStoreLeadCurrQCompPtr; nextStoreLeadCurrQCompPtr=0;}
  if (currStoreLeadCurrQCompDerivPtr) { delete currStoreLeadCurrQCompDerivPtr; currStoreLeadCurrQCompDerivPtr=0;}
  if (lastStoreLeadCurrQCompDerivPtr) { delete lastStoreLeadCurrQCompDerivPtr; lastStoreLeadCurrQCompDerivPtr=0;}
  //if (oldStoreLeadCurrQCompDerivPtr) { delete oldStoreLeadCurrQCompDerivPtr; oldStoreLeadCurrQCompDerivPtr=0;}
  if (nextStoreLeadCurrQCompDerivPtr) { delete nextStoreLeadCurrQCompDerivPtr; nextStoreLeadCurrQCompDerivPtr=0;}

  if (currSolutionDerivPtr) { delete currSolutionDerivPtr; currSolutionDerivPtr=0;}
  if (lastSolutionDerivPtr) { delete lastSolutionDerivPtr; lastSolutionDerivPtr=0;}

  if (nextSolutionDerivPtr) { delete nextSolutionDerivPtr; nextSolutionDerivPtr=0;}

  if (currStateDerivPtr) { delete currStateDerivPtr; currStateDerivPtr=0;}
  if (lastStateDerivPtr) { delete lastStateDerivPtr; lastStateDerivPtr=0;}

  if (nextStateDerivPtr) { delete nextStateDerivPtr;nextStateDerivPtr=0;}

  if (currSolutionDivDiffPtr) { delete currSolutionDivDiffPtr; currSolutionDivDiffPtr=0;}
  if (lastSolutionDivDiffPtr) { delete lastSolutionDivDiffPtr; lastSolutionDivDiffPtr=0;}

  if (nextSolutionDivDiffPtr) { delete nextSolutionDivDiffPtr; nextSolutionDivDiffPtr=0;}

  if (currStateDivDiffPtr) { delete currStateDivDiffPtr; currStateDivDiffPtr=0;}
  if (lastStateDivDiffPtr) { delete lastStateDivDiffPtr; lastStateDivDiffPtr=0;}

  if (nextStateDivDiffPtr) { delete nextStateDivDiffPtr; nextStateDivDiffPtr=0;}

  if (errWtVecPtr) { delete errWtVecPtr; errWtVecPtr=0;}
  if (absErrTolPtr) { delete absErrTolPtr; absErrTolPtr=0;}
  if (relErrTolPtr) { delete relErrTolPtr; relErrTolPtr=0;}

  if (deviceMaskPtr) { delete deviceMaskPtr; deviceMaskPtr=0;}

  if (newtonCorrectionPtr) { delete newtonCorrectionPtr; newtonCorrectionPtr=0;}

  //new-DAE:
   // Error Vectors
  if (qErrWtVecPtr) { delete qErrWtVecPtr; qErrWtVecPtr=0; }

  // DAE formulation vectors
  if (daeQVectorPtr) { delete daeQVectorPtr; daeQVectorPtr=0; }
  if (daeFVectorPtr) { delete daeFVectorPtr; daeFVectorPtr=0; }

  // DAE formulation matrices
  if (dQdxMatrixPtr) { delete dQdxMatrixPtr; dQdxMatrixPtr=0; }
  if (dFdxMatrixPtr) { delete dFdxMatrixPtr; dFdxMatrixPtr=0; }

  // HB temporary vectors
  if (dQdxVecVectorPtr) { delete dQdxVecVectorPtr; dQdxVecVectorPtr=0; }
  if (dFdxVecVectorPtr) { delete dFdxVecVectorPtr; dFdxVecVectorPtr=0; }

  // History arrays are STL containers of N_LAS_Vector 
  // so they'll take care of themselves.  NO!  they're containers of            
  // N_LAS_Vector *pointers*, so they need to be deleted by hand.               
  if (xHistory.size() != 0) 
  { 
    int sizeOfHistory=xHistory.size(); 
    for (int i=0;i<sizeOfHistory;++i) 
    { 
      if (xHistory[i]) {delete xHistory[i]; } 
      if (qHistory[i]) {delete qHistory[i]; } 
      if (sHistory[i]) {delete sHistory[i]; } 
      if (stoHistory[i]) {delete stoHistory[i]; }
      if (stoLeadCurrQCompHistory[i]) {delete stoLeadCurrQCompHistory[i];}
    } 
    xHistory.clear(); 
    qHistory.clear(); 
    sHistory.clear(); 
    stoHistory.clear();
    stoLeadCurrQCompHistory.clear();
  } 

  if (qn0Ptr)  { delete qn0Ptr; qn0Ptr=0; }
  if (qn0Ptr)  { delete qn0Ptr; qn0Ptr=0; }
  if (qpn0Ptr) { delete qpn0Ptr; qpn0Ptr=0; }
  if (sn0Ptr)  { delete sn0Ptr; sn0Ptr=0; }
  if (spn0Ptr) { delete spn0Ptr; spn0Ptr=0; }
  if (ston0Ptr)  { delete ston0Ptr; ston0Ptr=0; }
  if (stopn0Ptr) { delete stopn0Ptr; stopn0Ptr=0; }
  if (stoQCn0Ptr)  { delete stoQCn0Ptr; stoQCn0Ptr=0; }
  if (stoQCpn0Ptr) { delete stoQCpn0Ptr; stoQCpn0Ptr=0; }
  
  // Nonlinear solution vector:
  if (qNewtonCorrectionPtr) { delete qNewtonCorrectionPtr; qNewtonCorrectionPtr=0; }
  if (sNewtonCorrectionPtr) { delete sNewtonCorrectionPtr; sNewtonCorrectionPtr=0; }
  if (stoNewtonCorrectionPtr) { delete stoNewtonCorrectionPtr; stoNewtonCorrectionPtr=0; }
  if (stoLeadCurrQCompNewtonCorrectionPtr) { delete stoLeadCurrQCompNewtonCorrectionPtr; stoLeadCurrQCompNewtonCorrectionPtr=0; }
  
  // Step-size selection temporary vectors
  if (delta_x) { delete delta_x; delta_x=0; }
  if (delta_q) { delete delta_q; delta_q=0; }

  // Temporary vector for WaMPDE interpolation
  if (tmpXn0APtr) { delete tmpXn0APtr; tmpXn0APtr=0; }
  if (tmpXn0BPtr) { delete tmpXn0BPtr; tmpXn0BPtr=0; }

  // Delete data in the fast time storage for HB and MPDE
  resetFastTimeData();
  
  return ;
}

// ----------------------------------------------------------------------------
// -----------------------  DataStore Class Functions -------------------------
// ----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::printOutPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/23/00
//-----------------------------------------------------------------------------

void N_TIA_DataStore::printOutPointers()
{
  cout << "olde ptr = " << oldeSolutionPtr << endl;
  cout << "last ptr = " << lastSolutionPtr << endl;
  cout << "curr ptr = " << currSolutionPtr << endl;
  cout << "next ptr = " << nextSolutionPtr << endl;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setConstantHistory
// Purpose       : This function is called  after the operating point
//                 calculation has been called.  Once the operating point
//                 solution has been obtained, the code should regard that
//                 solution as having been the existing constant solution since
//                 the dawn of time.
// Special Notes : The most recent solution, etc., are in the "next" vectors.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/23/00
//-----------------------------------------------------------------------------
void N_TIA_DataStore::setConstantHistory()
{
#ifdef Xyce_DEBUG_TIME
  cout << "\nN_TIA_DataStore::setConstantHistory" << endl;
#endif

  // Solutions:
  *(lastSolutionPtr) = *(nextSolutionPtr);
  *(currSolutionPtr) = *(nextSolutionPtr);

  // Derivative of Solutions:
  *(lastSolutionDerivPtr) = *(nextSolutionDerivPtr);
  *(currSolutionDerivPtr) = *(nextSolutionDerivPtr);

  // Scaled Divided Differences of solutions:
  *(lastSolutionDivDiffPtr) = *(nextSolutionDivDiffPtr);
  *(currSolutionDivDiffPtr) = *(nextSolutionDivDiffPtr);

  // States:
  *(lastStatePtr) = *(nextStatePtr);
  *(currStatePtr) = *(nextStatePtr);

  // Stores:
  *(lastStorePtr) = *(nextStorePtr);
  *(currStorePtr) = *(nextStorePtr);
   
  *(lastStoreLeadCurrQCompPtr) = *(nextStoreLeadCurrQCompPtr);
  *(currStoreLeadCurrQCompPtr) = *(nextStoreLeadCurrQCompPtr);
  
  // Derivative of States:
  *(lastStateDerivPtr) = *(nextStateDerivPtr);
  *(currStateDerivPtr) = *(nextStateDerivPtr);
  
  // lead current derivative info in store 
  *(lastStoreLeadCurrQCompDerivPtr) = *(nextStoreLeadCurrQCompDerivPtr);
  *(currStoreLeadCurrQCompDerivPtr) = *(nextStoreLeadCurrQCompDerivPtr);
  
  // Scaled Divided Differences of states:
  *(lastStateDivDiffPtr) = *(nextStateDivDiffPtr);
  *(currStateDivDiffPtr) = *(nextStateDivDiffPtr);

}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::resetAll
//
// Purpose       : This function resets everything so that a transient loop
//                 can be started from the beginning.
//
// Special Notes : This function was needed for HB.
//
// Scope         : public
// Creator       : T. Mei, SNL 
// Creation Date : 02/26/09
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::resetAll ()
{
  absErrTolPtr->putScalar(tiaParamsPtr_->absErrorTol);
  relErrTolPtr->putScalar(tiaParamsPtr_->relErrorTol);
      
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::resetFastTimeData
//
// Purpose       : This function deletes the information from all the vectors
//                 that store fast time data for HB and MPDE
//
// Special Notes : This function was needed for HB.
//
// Scope         : public
// Creator       : Heidi Thornquist, SNL 
// Creation Date : 06/05/13
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::resetFastTimeData()
{
  // Clear the time step vectors
  timeSteps.clear();
  timeStepsBreakpointFlag.clear();

  // Delete any stored up any solution or state info
  std::vector<N_LAS_Vector *>::iterator currVecPtr = fastTimeSolutionVec.begin();
  std::vector<N_LAS_Vector *>::iterator endVecPtr = fastTimeSolutionVec.end();
  while( currVecPtr != endVecPtr )
  {
    delete *currVecPtr;
    currVecPtr++;
  }
  fastTimeSolutionVec.clear(); 
 
  currVecPtr = fastTimeStateVec.begin();
  endVecPtr = fastTimeStateVec.end();
  while( currVecPtr != endVecPtr )
  {
    delete *currVecPtr;
    currVecPtr++;
  }
  fastTimeStateVec.clear(); 
  
  currVecPtr = fastTimeQVec.begin();
  endVecPtr = fastTimeQVec.end();
  while( currVecPtr != endVecPtr )
  {
    delete *currVecPtr;
    currVecPtr++;
  }
  fastTimeQVec.clear(); 

  currVecPtr = fastTimeStoreVec.begin();
  endVecPtr = fastTimeStoreVec.end();
  while( currVecPtr != endVecPtr )
  {
    delete *currVecPtr;
    currVecPtr++;
  }
  fastTimeStoreVec.clear(); 
      
  return true;
}



//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::updateSolDataArrays
// Purpose       : Update the necessary integration data arrays for
//                 preparation of the next integration step. This is done after
//                 a successful step has been taken.
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
void N_TIA_DataStore::updateSolDataArrays()
{
#ifdef Xyce_DEBUG_TIME
  if (tiaParamsPtr_->debugLevel > 1)
  {
    cout << "\nN_TIA_DataStore::updateSolDataArrays " << endl;
  }
#endif

  // if the next solution has been switched out (probably because of NOX),
  // then reset it
  if (nextSolPtrSwitched) unsetNextSolVectorPtr();

  // Solutions:
  oldeSolutionPtr = lastSolutionPtr;
  lastSolutionPtr = currSolutionPtr;
  currSolutionPtr = nextSolutionPtr;
  nextSolutionPtr = oldeSolutionPtr;

  // Derivative of Solutions:
  oldeSolutionDerivPtr = lastSolutionDerivPtr;
  lastSolutionDerivPtr = currSolutionDerivPtr;
  currSolutionDerivPtr = nextSolutionDerivPtr;
  nextSolutionDerivPtr = oldeSolutionDerivPtr;

  // Scaled Divided Differences of solutions:
  oldeSolutionDivDiffPtr = lastSolutionDivDiffPtr;
  lastSolutionDivDiffPtr = currSolutionDivDiffPtr;
  currSolutionDivDiffPtr = nextSolutionDivDiffPtr;
  nextSolutionDivDiffPtr = oldeSolutionDivDiffPtr;

  // States:
  oldeStatePtr = lastStatePtr;
  lastStatePtr = currStatePtr;
  currStatePtr = nextStatePtr;
  nextStatePtr = oldeStatePtr;

  // Stores:
  oldeStorePtr = lastStorePtr;
  lastStorePtr = currStorePtr;
  currStorePtr = nextStorePtr;
  nextStorePtr = oldeStorePtr;
  oldStoreLeadCurrQCompPtr = lastStoreLeadCurrQCompPtr;
  lastStoreLeadCurrQCompPtr = currStoreLeadCurrQCompPtr;
  currStoreLeadCurrQCompPtr = nextStoreLeadCurrQCompPtr;
  nextStoreLeadCurrQCompPtr = oldStoreLeadCurrQCompPtr;

  // Derivative of States:
  oldeStateDerivPtr = lastStateDerivPtr;
  lastStateDerivPtr = currStateDerivPtr;
  currStateDerivPtr = nextStateDerivPtr;
  nextStateDerivPtr = oldeStateDerivPtr;
  
  // lead curent component of store 
  oldStoreLeadCurrQCompDerivPtr = lastStoreLeadCurrQCompDerivPtr;
  lastStoreLeadCurrQCompDerivPtr = currStoreLeadCurrQCompDerivPtr;
  currStoreLeadCurrQCompDerivPtr = nextStoreLeadCurrQCompDerivPtr;
  nextStoreLeadCurrQCompDerivPtr = oldStoreLeadCurrQCompDerivPtr;

  // Scaled Divided Differences of states:
  oldeStateDivDiffPtr = lastStateDivDiffPtr;
  lastStateDivDiffPtr = currStateDivDiffPtr;
  currStateDivDiffPtr = nextStateDivDiffPtr;
  nextStateDivDiffPtr = oldeStateDivDiffPtr;

  // copy contents of "curr" into "next".  This is to insure 
  // that at a minimum, the initial guess for the Newton solve 
  // will at least be the results of the previous Newton solve. 
  *(nextSolutionPtr) = *(currSolutionPtr);
  *(nextStatePtr)    = *(currStatePtr);
  *(nextStorePtr)    = *(currStorePtr);

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::updateStateDataArrays
//
//
// Purpose       : Same as updateSolDataArrays, but this function only 
//                 advances the state vector, and leaves the 
//                 solution alone.
//
// Special Notes : The main usage of this function is LOCA.  
//                 After each continuation step, LOCA needs to call 
//                 this function.  LOCA keeps track of solution vectors 
//                 on its own, which is why updateSolDataArrays
//                 is inappropriate for LOCA.
//
//                 This is necessary for voltage limiting to be 
//                 consistent with LOCA.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, 9233.
// Creation Date : 3/06/05
//-----------------------------------------------------------------------------

bool N_TIA_DataStore::updateStateDataArrays()
{
#ifdef Xyce_DEBUG_TIME
  cout << "\nN_TIA_DataStore::updateStateDataArrays " << endl;
#endif
  // States:
  oldeStatePtr = lastStatePtr;
  lastStatePtr = currStatePtr;
  currStatePtr = nextStatePtr;
  nextStatePtr = oldeStatePtr;

  // Stores:
  oldeStorePtr = lastStorePtr;
  lastStorePtr = currStorePtr;
  currStorePtr = nextStorePtr;
  nextStorePtr = oldeStorePtr;
  
  // q component of lead current 
  oldStoreLeadCurrQCompPtr = lastStoreLeadCurrQCompPtr;
  lastStoreLeadCurrQCompPtr = currStoreLeadCurrQCompPtr;
  currStoreLeadCurrQCompPtr = nextStoreLeadCurrQCompPtr;
  nextStoreLeadCurrQCompPtr = oldStoreLeadCurrQCompPtr;
  
  // Derivative of States:
  oldeStateDerivPtr = lastStateDerivPtr;
  lastStateDerivPtr = currStateDerivPtr;
  currStateDerivPtr = nextStateDerivPtr;
  nextStateDerivPtr = oldeStateDerivPtr;

  // Scaled Divided Differences of states:
  oldeStateDivDiffPtr = lastStateDivDiffPtr;
  lastStateDivDiffPtr = currStateDivDiffPtr;
  currStateDivDiffPtr = nextStateDivDiffPtr;
  nextStateDivDiffPtr = oldeStateDivDiffPtr;


  // Now, make the "next" stuff the same as the "curr" stuff.  
  // This is done because at the end of the tranop, but before
  // the transient phase starts, the function setConstantHistory
  // will be called.  When it is called, the most recent values
  // for state variables need to be in the "next" vectors.
  //
  // As long as LOCA solves are never used for transient, and only
  // for DC and tranop solves, this is OK.  This is, of course,
  // a bit of a kludge.  
  *(nextStatePtr) = *(currStatePtr);
  *(nextStorePtr) = *(currStorePtr);

  return true;
}



#ifdef Xyce_DEBUG_TIME
#ifndef OLD_PLOT
//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::outputSolDataArrays
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------

void N_TIA_DataStore::outputSolDataArrays()

{
  char tmp[256];
  static string dashedline =
"----------------------------------------------------------------------------";

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
 "  Solution Vectors:\n          Current               Last               Olde               Error");


  for (unsigned int k = 0; k < solutionSize; ++k)
  {
    sprintf(tmp,"%19.6e%19.6e%19.6e%19.6e",
      (*(currSolutionPtr))[k],
      (*(lastSolutionPtr))[k],
      (*(oldeSolutionPtr))[k],
      (*(newtonCorrectionPtr))[k]);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, tmp);
  }

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  return;
}

#else

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::outputSolDataArrays
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------

void N_TIA_DataStore::outputSolDataArrays()
{
  char tmp[128];

  for (unsigned int k = 0; k < solutionSize; ++k)
  {
    sprintf(tmp,"\t%14.6e", (*(currSolutionPtr))[0][k]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, tmp);
  }
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
}

#endif
#endif // Xyce_DEBUG_TIME

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::enableOrderOneStart
// Purpose       : Initialize arrays for "re-starting" integration at order 1.
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------

void N_TIA_DataStore::enableOrderOneStart()
{
  // solution vectors:
  *(lastSolutionPtr) = *(currSolutionPtr);
  *(oldeSolutionPtr) = *(currSolutionPtr);

  nextSolutionDerivPtr->putScalar(0.0);
  currSolutionDerivPtr->putScalar(0.0);
  currSolutionDivDiffPtr->putScalar(0.0);

  // state vectors:
  *(lastStatePtr) = *(currStatePtr);
  *(oldeStatePtr) = *(currStatePtr);

  nextStateDerivPtr->putScalar(0.0);
  currStateDerivPtr->putScalar(0.0);
  currStateDivDiffPtr->putScalar(0.0);

  // store vectors:
  *(lastStorePtr) = *(currStorePtr);
  *(oldeStorePtr) = *(currStorePtr);
  *(lastStoreLeadCurrQCompPtr) = *(currStoreLeadCurrQCompPtr);
  *(oldStoreLeadCurrQCompPtr)= *(currStoreLeadCurrQCompPtr);
  nextStoreLeadCurrQCompDerivPtr->putScalar(0.0);
  currStoreLeadCurrQCompDerivPtr->putScalar(0.0);
  
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::outputPredictedSolution
// Purpose       : Set ErrorEstimate array values to zero.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
#ifdef Xyce_DEBUG_TIME
void N_TIA_DataStore::outputPredictedSolution()

{
  char tmp[128];
  static string dashedline2 = "---------------------";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,dashedline2);

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,"  Predicted Solution:");

  for (unsigned int k = 0; k< solutionSize; ++k)
  {
    sprintf(tmp,"%14.6e", (*(nextSolutionPtr))[k]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,tmp);
  }

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,dashedline2);
}
#endif

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::outputPredictedDerivative
// Purpose       : Set ErrorEstimate array values to zero.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
#ifdef Xyce_DEBUG_TIME
void N_TIA_DataStore::outputPredictedDerivative()
{
  char tmp[128];
  static string dashedline2 = "---------------------";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,dashedline2);

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,"  Predicted Derivative:");


  for (unsigned int k = 0; k < solutionSize; ++k)
  {
    sprintf(tmp,"%14.6e", (*(nextSolutionDerivPtr))[k] );
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,tmp);
  }

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,dashedline2);
}
#endif

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::partialErrorNormSum
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
// 
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2 
//
//                  It will later be summed with other patial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
double N_TIA_DataStore::partialErrorNormSum ()
{
  double errorNorm = 0.0;
  newtonCorrectionPtr->wRMSNorm(*errWtVecPtr, &errorNorm);
  double sum = errorNorm*errorNorm;
  double length = newtonCorrectionPtr->globalLength();

  sum *= length;

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::globalLength
// Purpose       : Needed by 2-level solves.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
double N_TIA_DataStore::globalLength ()
{
  return newtonCorrectionPtr->globalLength();
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::computeDividedDifferences
// Purpose       : Compute the scaled (by current stepsize) divided difference
//                 approximation to the derivative.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
void N_TIA_DataStore::computeDividedDifferences()
{
  // First do the solution divided difference:
  nextSolutionDivDiffPtr->linearCombo(1.0, *(nextSolutionPtr), -1.0, *(currSolutionPtr));

  // Now do the state divided difference:
  nextStateDivDiffPtr->linearCombo(1.0, *(nextStatePtr), -1.0, *(currStatePtr));

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::computeDivDiffsBlock
// Purpose       : Compute the scaled (by current stepsize) divided difference
//                 approximation to the derivative.  This is the same as
//                 the  function "computeDividedDifferences", except that
//                 the  operation is only performed on a sub-block of the
//                 vectors.
//
// Special Notes : Not done yet...
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/10/01
//-----------------------------------------------------------------------------
void N_TIA_DataStore::computeDivDiffsBlock (
         const list<index_pair> & solGIDList,
         const list<index_pair> & staGIDList)
{

#if 0
  computeDividedDifferences();
#endif

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::equateTmpVectors
// Purpose       : This function equates the 6 temporary vectors  with
//                 their "next" vector equivalents.  This function
//                 is neccessary for the nonlinear solver damping loop.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/02
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::equateTmpVectors()
{
  // next solution vector:
  *(tmpSolVectorPtr) = *(nextSolutionPtr);

  // next state vector:
  *(tmpStaVectorPtr) = *(nextStatePtr);

  // next store vector:
  *(tmpStoVectorPtr) = *(nextStorePtr);

  // next state divided difference vector:
  *(tmpStaDivDiffPtr) = *(nextStateDivDiffPtr);

  // next state derivative vector:
  *(tmpStaDerivPtr)  = *(nextStateDerivPtr);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::usePreviousSolAsPredictor
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/28/01
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::usePreviousSolAsPredictor ()
{
  bool bsuccess = true;

  *(nextSolutionPtr) = *(currSolutionPtr);
  *(nextStatePtr)    = *(currStatePtr);
  *(nextStorePtr)    = *(currStorePtr);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setNextSolVectorPtr
// Purpose       :
// Special Notes : Only needed for NOX, and other solvers that prefer to
//                 own the solution vector.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/03
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::setNextSolVectorPtr (N_LAS_Vector * solVecPtr)
{
  // only save the old pointer if it hasn't been switched yet.
  if (!nextSolPtrSwitched)
  {
    savedNextSolutionPtr = nextSolutionPtr;
    nextSolPtrSwitched = true;
  }
  nextSolutionPtr = solVecPtr;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::unsetNextSolVectorPtr
// Purpose       :
// Special Notes : This also copies over the solution, in addition to the
//                 pointer.
//
//                 This is only called when it is time to rotate the
//                 pointers for the next time step.  If we have been
//                 running with NOX, or with some other solver that prefers
//                 to own the solution vector, this is neccessary, and the
//                 switch flag should be set.  Otherwise, not.
//
//                 Basically, if we've been running with NOX, then the next
//                 solution vector ptr has probably been switched out at least
//                 once.  We need to maintain the history, so we make a
//                 copy of this switched solution vector, and the 
//                 restore the old pointer.
//
//                 This is kludgy, but will have to do for now.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/03
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::unsetNextSolVectorPtr ()
{
  if (nextSolPtrSwitched)
  {
    *(savedNextSolutionPtr) = *(nextSolutionPtr);
    nextSolutionPtr = savedNextSolutionPtr;
    nextSolPtrSwitched = false;
  }
  return true;
}

// new-DAE data store functions:


//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setZeroHistory
// Purpose       : Sets everything to zero.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/24/07
//-----------------------------------------------------------------------------
void N_TIA_DataStore::setZeroHistory()
{
#ifdef Xyce_DEBUG_TIME
  cout << "\nN_TIA_DataStore::setZeroHistory" << endl;
#endif
  // Solutions:
  nextSolutionPtr->putScalar(0.0);
  nextSolutionDerivPtr->putScalar(0.0);
  nextSolutionDivDiffPtr->putScalar(0.0);

  // States:
  nextStatePtr->putScalar(0.0);
  nextStateDerivPtr->putScalar(0.0);
  nextStateDivDiffPtr->putScalar(0.0);
  nextStorePtr->putScalar(0.0);
  nextStoreLeadCurrQCompPtr->putScalar(0.0);

  qErrWtVecPtr->putScalar(0.0);

  // DAE formulation vectors
  daeQVectorPtr->putScalar(0.0);
  daeFVectorPtr->putScalar(0.0);

  // Predictors
  xn0Ptr->putScalar(0.0);
  qn0Ptr->putScalar(0.0);
  qpn0Ptr->putScalar(0.0);
  sn0Ptr->putScalar(0.0);
  spn0Ptr->putScalar(0.0);
  stoQCn0Ptr->putScalar(0.0);
  stoQCpn0Ptr->putScalar(0.0);

  // Nonlinear solution vector:
  qNewtonCorrectionPtr->putScalar(0.0);
  sNewtonCorrectionPtr->putScalar(0.0);
  stoNewtonCorrectionPtr->putScalar(0.0);
  stoLeadCurrQCompNewtonCorrectionPtr->putScalar(0.0);

  // Step-size selection temporary vectors
  delta_x->putScalar(0.0);
  delta_q->putScalar(0.0);

  // Temporary vector for WaMPDE interpolation
  tmpXn0APtr->putScalar(0.0);
  tmpXn0BPtr->putScalar(0.0);

  // This just sets the "oldDAE" history vectors to zero.
  setConstantHistory ();

  // new-DAE history:
  int sizeOfHistory = xHistory.size();
  for (int i=0;i<sizeOfHistory;++i)
  {
    (xHistory[i])->putScalar(0.0);
    (qHistory[i])->putScalar(0.0);
    (sHistory[i])->putScalar(0.0);
    (stoHistory[i])->putScalar(0.0);
    (stoLeadCurrQCompHistory[i])->putScalar(0.0);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setErrorWtVector_
// Purpose       : Set the Error Weight Vector (defined in terms of the
//                 solution approximation and error tolerances).
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL,Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
void N_TIA_DataStore::setErrorWtVector()
{
  // to avoid several conditionals within a loop we traverse many times, we'll
  // presort the variables into types so we can loop without also making conditional
  // checks
  if( ! indexVecsInitialized )
  {
    // figure out which unknowns are V's, I's or masked.
    numVVars = 0;
    numIVars = 0;
    numMaskedVars = 0;
    
    bool nTDMF=lasSysPtr->getNonTrivialDeviceMaskFlag();
    
    // first count how many of each type we have
    for (int k = 0; k < solutionSize; ++k)
    {
      if (nTDMF && (*(deviceMaskPtr))[k] == 0.0)
      {
        numMaskedVars++;
      }
      else if( (tiaParamsPtr_->fastTests == true) && (varTypeVec[k] == 'I') )
      {
        // we only count I vars if we're doing fast tests, otherwise we treat them as v vars
        numIVars++;
      }
      else
      {
        numVVars++;   
      }
    }
    
    // set up our storage
    indexVVars.resize(numVVars);
    indexIVars.resize(numIVars);
    indexMaskedVars.resize(numMaskedVars);

    // now fill the index arrays
    for (int k = 0, currI = 0, currV = 0, currM = 0 ; k < solutionSize; ++k)
    {
      if (nTDMF && (*(deviceMaskPtr))[k] == 0.0)
      {
        indexMaskedVars[currM] = k;
        currM++;
      }
      else if( (tiaParamsPtr_->fastTests == true) && (varTypeVec[k] == 'I') )
      {
        indexIVars[currI] = k;
        currI++;
      }
      else
      {
        indexVVars[currV] = k;
        currV++;
      }
    }
    
    indexVecsInitialized = true;
  }

#ifdef Xyce_DEBUG_TIME
  static string dashedline =
"----------------------------------------------------------------------------";
  char tmp[256];

  if (tiaParamsPtr_->debugLevel > 1)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "N_TIA_DataStore::setErrorWtVector");

    sprintf(tmp,
    "\n   errorWtVector    currSolution     relErrorTol    absErrorTol");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, tmp);

    sprintf(tmp,
    "   --------------  --------------  --------------  --------------");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, tmp);
  }
#endif

  if (tiaParamsPtr_->newLte == true)
//  if (tiaParamsPtr_->integrationMethod == 7 && tiaParamsPtr_->newLte == true)
  { 
    double currMaxValue = 0.0;
    currSolutionPtr->infNorm(&currMaxValue);
 
    errWtVecPtr->putScalar(currMaxValue);

#ifdef Xyce_DEBUG_TIME 
    if (tiaParamsPtr_->debugLevel > 0)
      {
        std::vector<int> index(1, -1);
        currSolutionPtr->infNormIndex( &index[0] );
        std::string outMessage = "currMaxValue = " + Teuchos::Utils::toString( currMaxValue )
                                 + ", currMaxValueIndex = " + Teuchos::Utils::toString( index[0] );
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, outMessage);
      }
#endif // Xyce_DEBUG_TIME 
  }
  else
  {
    errWtVecPtr->absValue(*currSolutionPtr);

#ifdef Xyce_DEBUG_TIME 
    if (tiaParamsPtr_->debugLevel > 0)
      {
        double currMaxValue = 0.0;
        currSolutionPtr->infNorm(&currMaxValue);
        std::vector<int> index(1, -1);
        currSolutionPtr->infNormIndex( &index[0] );
        std::string outMessage = "currMaxValueoldLte = " + Teuchos::Utils::toString( currMaxValue )
                                 + ", currMaxValueIndex = " + Teuchos::Utils::toString( index[0] );
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, outMessage);
      }
#endif
    
  }

  qErrWtVecPtr->absValue(*daeQVectorPtr);

  if( tiaParamsPtr_->fastTests == true )
  {
    // Voltage variables
    for (int k = 0; k < numVVars; ++k)
    {
//      errWtVecPtr->putScalar(currMaxValue);
      if( (*(errWtVecPtr))[indexVVars[k]] < tiaParamsPtr_->voltZeroTol )
      { 
        (*(errWtVecPtr))[indexVVars[k]] = N_UTL_MachineDependentParams::MachineBig();
        (*(qErrWtVecPtr))[indexVVars[k]] = N_UTL_MachineDependentParams::MachineBig();
      }
      else
      {
        (*(errWtVecPtr))[indexVVars[k]] = (*(relErrTolPtr))[indexVVars[k]] * (*(errWtVecPtr))[indexVVars[k]] + (*(absErrTolPtr))[indexVVars[k]];  
        (*(qErrWtVecPtr))[indexVVars[k]] = (*(relErrTolPtr))[indexVVars[k]] * (*(qErrWtVecPtr))[indexVVars[k]] + (*(absErrTolPtr))[indexVVars[k]];
      }
    }
    
    
    // Current variables
    for (int k = 0; k < numIVars; ++k)
    {
//      errWtVecPtr->absValue(*currSolutionPtr);
      if( (*(errWtVecPtr))[indexIVars[k]] < tiaParamsPtr_->currZeroTol )
      { 
        (*(errWtVecPtr))[indexIVars[k]] = N_UTL_MachineDependentParams::MachineBig();
        (*(qErrWtVecPtr))[indexIVars[k]] = N_UTL_MachineDependentParams::MachineBig();
      }
      else
      {
        (*(errWtVecPtr))[indexIVars[k]] =  (*(errWtVecPtr))[indexIVars[k]] + (*(absErrTolPtr))[indexIVars[k]];
        (*(qErrWtVecPtr))[indexIVars[k]] = (*(qErrWtVecPtr))[indexIVars[k]] + (*(absErrTolPtr))[indexIVars[k]];
      }
    }
    
    // Masked variables
    for (int k = 0; k < numMaskedVars; ++k)
    {
      (*(errWtVecPtr))[indexMaskedVars[k]] = (*(qErrWtVecPtr))[indexMaskedVars[k]] = N_UTL_MachineDependentParams::MachineBig();
    }
  }
  else
  {
     // Voltage variables
    for (int k = 0; k < numVVars; ++k)
    {
//      errWtVecPtr->putScalar(currMaxValue);
      (*(errWtVecPtr))[indexVVars[k]] = (*(relErrTolPtr))[indexVVars[k]] * (*(errWtVecPtr))[indexVVars[k]] + (*(absErrTolPtr))[indexVVars[k]];  
      (*(qErrWtVecPtr))[indexVVars[k]] = (*(relErrTolPtr))[indexVVars[k]] * (*(qErrWtVecPtr))[indexVVars[k]] + (*(absErrTolPtr))[indexVVars[k]];
    }
    
    // Current variables
    // if fastTests == false, then I vars are treated with V vars above.
    
    // Masked variables
    for (int k = 0; k < numMaskedVars; ++k)
    {
      (*(errWtVecPtr))[indexMaskedVars[k]] = (*(qErrWtVecPtr))[indexMaskedVars[k]] = N_UTL_MachineDependentParams::MachineBig();
    }
 
  }

#ifdef Xyce_DEBUG_TIME
  for (int k = 0; k < solutionSize; ++k)
  {
    if (tiaParamsPtr_->debugLevel > 1)
    {
      sprintf(tmp,"%16.6e%16.6e%16.6e%16.6e",
              (*(errWtVecPtr))     [k],
              (*(currSolutionPtr)) [k],
              (*(relErrTolPtr))    [k],
              (*(absErrTolPtr))    [k]);

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, tmp);
    }
  }
#endif

#ifdef Xyce_DEBUG_TIME
  if (tiaParamsPtr_->debugLevel > 1)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "" );
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::WRMS_errorNorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
double N_TIA_DataStore::WRMS_errorNorm()
{
  double errorNorm = 0.0, qErrorNorm = 0.0;
  newtonCorrectionPtr->wRMSNorm(*errWtVecPtr, &errorNorm);
  qNewtonCorrectionPtr->wRMSNorm(*qErrWtVecPtr, &qErrorNorm);

#ifdef Xyce_EXTDEV
#ifdef Xyce_DEBUG_TIME
  if (tiaParamsPtr_->debugLevel > 1)
  {
    cout << "N_TIA_DataStore::errorNorm = " << errorNorm << endl;
    cout << "N_TIA_DataStore::qErrorNorm = " << qErrorNorm << endl;
  }
#endif

  // This is for handling the total errorNorm for 2-level solves.
  //
  // Note:  This version of the function assumes that the error norm
  // and q-error norm are the same size.  (they have to be...)
  if ( !(innerErrorInfoVec.empty()) )
  {
    double upperSize = newtonCorrectionPtr->globalLength();
    int sumSize = innerErrorInfoVec.size();

    double totalSize = upperSize;
    double totalSum = errorNorm*errorNorm*upperSize;
    double totalQSum = qErrorNorm*qErrorNorm*upperSize;

    for (int i=0;i<sumSize;++i)
    {
      double innerSum = innerErrorInfoVec[i].xErrorSum;
      double innerQSum = innerErrorInfoVec[i].qErrorSum;
      double innerSize = innerErrorInfoVec[i].innerSize;

#ifdef Xyce_DEBUG_TIME
      cout << "DSdae:innerSum["<<i<<"] = " << innerSum <<endl;
      cout << "DSdae:innerQSum["<<i<<"] = " << innerQSum <<endl;
      cout << "DSdae:innerSize["<<i<<"] = " << innerSize <<endl;
#endif

      totalSize += innerSize;
      totalSum += innerSum;
      totalQSum += innerQSum;
    }

    double recip = 1.0/totalSize;
    errorNorm = sqrt(recip*totalSum);
    qErrorNorm = sqrt(recip*totalQSum);

#ifdef Xyce_DEBUG_TIME
    cout << "DSdae:upperSize = " << upperSize << endl;
    cout << "DSdae:totalSum = " << totalSum << endl;
    cout << "DSdae:totalQSum = " << totalQSum << endl;
    cout << "DSdae:totalSize = " << totalSize << endl;
    cout << "DSdae:2-level errorNorm = " << errorNorm << endl;
    cout << "DSdae:2-level qErrorNorm = " << qErrorNorm << endl;
#endif
  }
#endif

#ifndef Xyce_USE_Q_NORM
  //errorNorm = errorNorm;
#else
  errorNorm = sqrt(0.5*errorNorm*errorNorm+0.5*qErrorNorm*qErrorNorm);
#endif
  return errorNorm;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::partialQErrorNormSum
// Purpose       : Needed by 2-level solves.  This is the Q-vector version
//                 of this function.
//
// Special Notes : A weighted RMS norm is this:
// 
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2 
//
//                  It will later be summed with other patial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
double N_TIA_DataStore::partialQErrorNormSum ()
{
  double qErrorNorm = 0.0;
  qNewtonCorrectionPtr->wRMSNorm(*qErrWtVecPtr, &qErrorNorm);
  double sum = qErrorNorm*qErrorNorm;
  double length = qNewtonCorrectionPtr->globalLength();
  sum *= length;

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::partialSum_m1
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
// 
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2 
//
//                  It will later be summed with other patial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::partialSum_m1 (int currentOrder)
{
  double sum = 0.0;

  if (currentOrder>1)
  {
    delta_x->linearCombo(1.0,*(xHistory[currentOrder]),1.0,*newtonCorrectionPtr);
    double norm = 0.0;
    delta_x->wRMSNorm(*errWtVecPtr, &norm);
    sum = norm*norm;
    double length = newtonCorrectionPtr->globalLength();
    sum *= length;
  }

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::partialSum_m2
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
// 
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2 
//
//                  It will later be summed with other patial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::partialSum_m2 (int currentOrder)
{
  double sum = 0.0;

  if (currentOrder>2)
  {
    delta_x->linearCombo(1.0,*(xHistory[currentOrder]),1.0,*newtonCorrectionPtr);
    delta_x->linearCombo(1.0,*(xHistory[currentOrder-1]),1.0,*delta_x);
    double norm = 0.0;
    delta_x->wRMSNorm(*errWtVecPtr, &norm);
    sum = norm*norm;
    double length = newtonCorrectionPtr->globalLength();
    sum *= length;
  }

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::partialSum_p1
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
// 
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2 
//
//                  It will later be summed with other patial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::partialSum_p1 (int currentOrder, int maxOrder)
{
  double sum = 0.0;

  if (currentOrder<maxOrder)
  {
    delta_x->linearCombo(1.0,*newtonCorrectionPtr,-1.0,*(xHistory[currentOrder+1]));
    double norm = 0.0;
    delta_x->wRMSNorm(*errWtVecPtr, &norm);
    sum = norm*norm;
    double length = newtonCorrectionPtr->globalLength();
    sum *= length;
  }

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::partialSum_q1
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
// 
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2 
//
//                  It will later be summed with other patial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::partialSum_q1 ()
{
  double sum = 0.0;

  double norm = 0.0;
  (qHistory[1])->wRMSNorm(*qErrWtVecPtr, &norm);

  sum = norm*norm;
  double length = qNewtonCorrectionPtr->globalLength();
  sum *= length;

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::delta_x_errorNorm_m1
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::delta_x_errorNorm_m1()
{
  double errorNorm = 0.0;
  delta_x->wRMSNorm(*errWtVecPtr, &errorNorm);

#ifdef Xyce_EXTDEV
  // This is for handling the total errorNorm for 2-level solves.
  // It really should get a separate overloaded function.
  if ( !(innerErrorInfoVec.empty()) )
  {
    double upperSize = newtonCorrectionPtr->globalLength();
    int sumSize = innerErrorInfoVec.size();

    double totalSize = upperSize;
    double totalSum = errorNorm*errorNorm*upperSize;

    for (int i=0;i<sumSize;++i)
    {
      double innerSum = innerErrorInfoVec[i].xErrorSum_m1;
      double innerSize = innerErrorInfoVec[i].innerSize;

      totalSize += innerSize;
      totalSum += innerSum;
    }

    double recip = 1.0/totalSize;
    errorNorm = sqrt(recip*totalSum);
  }
#endif

  return errorNorm;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::delta_x_errorNorm_m2
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::delta_x_errorNorm_m2()
{
  double errorNorm = 0.0;
  delta_x->wRMSNorm(*errWtVecPtr, &errorNorm);

#ifdef Xyce_EXTDEV
  // This is for handling the total errorNorm for 2-level solves.
  // It really should get a separate overloaded function.
  if ( !(innerErrorInfoVec.empty()) )
  {
    double upperSize = newtonCorrectionPtr->globalLength();
    int sumSize = innerErrorInfoVec.size();

    double totalSize = upperSize;
    double totalSum = errorNorm*errorNorm*upperSize;

    for (int i=0;i<sumSize;++i)
    {
      double innerSum = innerErrorInfoVec[i].xErrorSum_m2;
      double innerSize = innerErrorInfoVec[i].innerSize;

      totalSize += innerSize;
      totalSum += innerSum;
    }

    double recip = 1.0/totalSize;
    errorNorm = sqrt(recip*totalSum);
  }
#endif

  return errorNorm;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::delta_x_errorNorm_p1
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::delta_x_errorNorm_p1()
{
  double errorNorm = 0.0;
  delta_x->wRMSNorm(*errWtVecPtr, &errorNorm);

#ifdef Xyce_EXTDEV
  // This is for handling the total errorNorm for 2-level solves.
  // It really should get a separate overloaded function.
  if ( !(innerErrorInfoVec.empty()) )
  {
    double upperSize = newtonCorrectionPtr->globalLength();
    int sumSize = innerErrorInfoVec.size();

    double totalSize = upperSize;
    double totalSum = errorNorm*errorNorm*upperSize;

    for (int i=0;i<sumSize;++i)
    {
      double innerSum = innerErrorInfoVec[i].xErrorSum_p1;
      double innerSize = innerErrorInfoVec[i].innerSize;

      totalSize += innerSize;
      totalSum += innerSum;
    }

    double recip = 1.0/totalSize;
    errorNorm = sqrt(recip*totalSum);
  }
#endif

  return errorNorm;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::delta_x_errorNorm_q1
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::delta_x_errorNorm_q1()
{
  double errorNorm = 0.0;
  (qHistory[1])->wRMSNorm(*qErrWtVecPtr, &errorNorm);

#ifdef Xyce_EXTDEV

  // This is for handling the total errorNorm for 2-level solves.
  // It really should get a separate overloaded function.
  if ( !(innerErrorInfoVec.empty()) )
  {
    double upperSize = newtonCorrectionPtr->globalLength();
    int sumSize = innerErrorInfoVec.size();

    double totalSize = upperSize;
    double totalSum = errorNorm*errorNorm*upperSize;

    for (int i=0;i<sumSize;++i)
    {
      double innerSum = innerErrorInfoVec[i].q1HistorySum;
      double innerSize = innerErrorInfoVec[i].innerSize;

      totalSize += innerSize;
      totalSum += innerSum;
    }

    double recip = 1.0/totalSize;
    errorNorm = sqrt(recip*totalSum);
  }

#endif

  return errorNorm;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::stepLinearCombo
// Purpose       : setup the newtonCorrection vectors.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 02/20/07
//-----------------------------------------------------------------------------
void N_TIA_DataStore::stepLinearCombo()
{
  // 03/16/04 tscoffe:  update the newton correction.  Note:  this should be
  // available from NOX, but for now I'm going to do the difference anyway.
  newtonCorrectionPtr->linearCombo (1.0,*nextSolutionPtr,-1.0,*xn0Ptr);

  // We need to compute the correction in Q here 
  // I'm assuming dsDaePtr_->daeQVectorPtr will be fresh from the end of the
  // nonlinear solve.
  qNewtonCorrectionPtr->linearCombo (1.0,*daeQVectorPtr,-1.0,*qn0Ptr);

  // We also need a State correction between the time steps
  sNewtonCorrectionPtr->linearCombo (1.0,*nextStatePtr,-1.0,*sn0Ptr);
  
  // We also need a Store correction between the time steps
  stoNewtonCorrectionPtr->linearCombo (1.0,*nextStorePtr,-1.0,*ston0Ptr);
  stoLeadCurrQCompNewtonCorrectionPtr->linearCombo (1.0, *nextStoreLeadCurrQCompPtr, -1.0, *stoQCn0Ptr);

#ifdef Xyce_DEBUG_TIME
  if (tiaParamsPtr_->debugLevel > 1)
  {
    cout << "\n newtonCorrection: \n" << endl;
    newtonCorrectionPtr->printPetraObject();
    cout << endl;
    cout << "\n qNewtonCorrection: \n" << endl;
    qNewtonCorrectionPtr->printPetraObject();
    cout << "\n sNewtonCorrection: \n" << endl;
    sNewtonCorrectionPtr->printPetraObject();
    cout << endl;
  }
#endif // Xyce_DEBUG_TIME

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::getSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::getSolnVarData( const int & gid,
				      vector<double> & varData )
{
  varData.resize(22);
  int i=0;
  varData[i++] = tmpSolVectorPtr->getElementByGlobalIndex( gid );
  varData[i++] = currSolutionPtr->getElementByGlobalIndex( gid );
  varData[i++] = lastSolutionPtr->getElementByGlobalIndex( gid );
  varData[i++] = nextSolutionPtr->getElementByGlobalIndex( gid );
  varData[i++] = currSolutionDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = lastSolutionDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = nextSolutionDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = currSolutionDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = lastSolutionDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = nextSolutionDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = errWtVecPtr->getElementByGlobalIndex( gid );
  varData[i++] = absErrTolPtr->getElementByGlobalIndex( gid );
  varData[i++] = relErrTolPtr->getElementByGlobalIndex( gid );
  varData[i++] = newtonCorrectionPtr->getElementByGlobalIndex( gid );
  varData[i++] = qErrWtVecPtr->getElementByGlobalIndex ( gid );
  varData[i++] = daeQVectorPtr->getElementByGlobalIndex ( gid );
  varData[i++] = daeFVectorPtr->getElementByGlobalIndex ( gid );
  varData[i++] = xn0Ptr->getElementByGlobalIndex ( gid );
  varData[i++] = qn0Ptr->getElementByGlobalIndex ( gid );
  varData[i++] = qpn0Ptr->getElementByGlobalIndex ( gid );
  varData[i++] = dFdxdVpVectorPtr->getElementByGlobalIndex ( gid );
  varData[i++] = dQdxdVpVectorPtr->getElementByGlobalIndex ( gid );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::getStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::getStateVarData( const int & gid,
			 	       vector<double> & varData )
{
  int i=0;
  varData.resize( 14 );
  varData[i++] = tmpStaVectorPtr->getElementByGlobalIndex( gid );
  varData[i++] = tmpStaDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = tmpStaDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = currStatePtr->getElementByGlobalIndex( gid );
  varData[i++] = lastStatePtr->getElementByGlobalIndex( gid );
  varData[i++] = nextStatePtr->getElementByGlobalIndex( gid );
  varData[i++] = currStateDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = lastStateDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = nextStateDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = currStateDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = lastStateDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = nextStateDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = sn0Ptr->getElementByGlobalIndex( gid );
  varData[i++] = spn0Ptr->getElementByGlobalIndex( gid );
  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::getStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::getStoreVarData( const int & gid,
			 	       vector<double> & varData )
{
  int i=0;
  varData.resize( 6 );
  varData[i++] = tmpStoVectorPtr->getElementByGlobalIndex( gid );
  varData[i++] = currStorePtr->getElementByGlobalIndex( gid );
  varData[i++] = lastStorePtr->getElementByGlobalIndex( gid );
  varData[i++] = nextStorePtr->getElementByGlobalIndex( gid );
  varData[i++] = ston0Ptr->getElementByGlobalIndex( gid );
  varData[i++] = stopn0Ptr->getElementByGlobalIndex( gid );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::setSolnVarData( const int & gid,
				      const vector<double> & varData )
{
  int i=0;
  tmpSolVectorPtr->setElementByGlobalIndex       ( gid, varData[i++] );
  currSolutionPtr->setElementByGlobalIndex       ( gid, varData[i++] );
  lastSolutionPtr->setElementByGlobalIndex       ( gid, varData[i++] );
  nextSolutionPtr->setElementByGlobalIndex       ( gid, varData[i++] );
  currSolutionDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );
  lastSolutionDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );
  nextSolutionDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );
  currSolutionDivDiffPtr->setElementByGlobalIndex( gid, varData[i++] );
  lastSolutionDivDiffPtr->setElementByGlobalIndex( gid, varData[i++] );
  nextSolutionDivDiffPtr->setElementByGlobalIndex( gid, varData[i++] );
  errWtVecPtr->setElementByGlobalIndex           ( gid, varData[i++] );
  absErrTolPtr->setElementByGlobalIndex          ( gid, varData[i++] );
  relErrTolPtr->setElementByGlobalIndex          ( gid, varData[i++] );
  newtonCorrectionPtr->setElementByGlobalIndex   ( gid, varData[i++] );
  qErrWtVecPtr->setElementByGlobalIndex          ( gid, varData[i++] );
  daeQVectorPtr->setElementByGlobalIndex         ( gid, varData[i++] );
  daeFVectorPtr->setElementByGlobalIndex         ( gid, varData[i++] );
  xn0Ptr->setElementByGlobalIndex                ( gid, varData[i++] );
  qn0Ptr->setElementByGlobalIndex                ( gid, varData[i++] );
  qpn0Ptr->setElementByGlobalIndex               ( gid, varData[i++] );
  dFdxdVpVectorPtr->setElementByGlobalIndex      ( gid, varData[i++] );
  dQdxdVpVectorPtr->setElementByGlobalIndex      ( gid, varData[i++] );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::setStateVarData( const int & gid,
			 	       const vector<double> & varData )
{
  int i=0;
  tmpStaVectorPtr->setElementByGlobalIndex    ( gid, varData[i++] );
  tmpStaDerivPtr->setElementByGlobalIndex     ( gid, varData[i++] );
  tmpStaDivDiffPtr->setElementByGlobalIndex   ( gid, varData[i++] );
  currStatePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  lastStatePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  nextStatePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  currStateDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );
  lastStateDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );
  nextStateDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );
  currStateDivDiffPtr->setElementByGlobalIndex( gid, varData[i++] );
  lastStateDivDiffPtr->setElementByGlobalIndex( gid, varData[i++] );
  nextStateDivDiffPtr->setElementByGlobalIndex( gid, varData[i++] );

  sn0Ptr->setElementByGlobalIndex             ( gid, varData[i++] );
  spn0Ptr->setElementByGlobalIndex            ( gid, varData[i++] );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::setStoreVarData( const int & gid,
			 	       const vector<double> & varData )
{
  int i=0;
  tmpStoVectorPtr->setElementByGlobalIndex    ( gid, varData[i++] );
  currStorePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  lastStorePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  nextStorePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  ston0Ptr->setElementByGlobalIndex           ( gid, varData[i++] );
  stopn0Ptr->setElementByGlobalIndex          ( gid, varData[i++] );
  return true;
}

