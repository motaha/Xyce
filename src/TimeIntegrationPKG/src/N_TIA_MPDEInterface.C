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
// Filename      : $RCSfile: N_TIA_MPDEInterface.C,v $
// Purpose       : This file contains the functions which define the time
//                 integration interface for the MPDE classes.
// Special Notes :
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/25/04
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.25.6.2 $
// Revision Date  : $Date: 2013/10/03 17:23:49 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <iostream>

// ----------   Xyce Includes   ----------

#include <N_TIA_MPDEInterface.h>
#include <N_TIA_TimeIntInfo.h> 
#include <N_TIA_TwoLevelError.h>

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::N_TIA_MPDEInterface
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/25/04
//-----------------------------------------------------------------------------
N_TIA_MPDEInterface::N_TIA_MPDEInterface(N_TIA_TIAParams & tp)
    : anaManagerPtr_(NULL),
      tiaParams_(tp),
      dsPtr_(NULL)
{
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::~N_TIA_MPDEInterface
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/25/04
//-----------------------------------------------------------------------------
N_TIA_MPDEInterface::~N_TIA_MPDEInterface ()
{
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::registerTIAParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/25/04
//-----------------------------------------------------------------------------
bool N_TIA_MPDEInterface::registerTIAParams(const N_TIA_TIAParams & tiaParams_tmp)

{
  tiaParams_ = tiaParams_tmp;
  if (anaManagerPtr_)
  {
    anaManagerPtr_->registerTIAParams (tiaParams_);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::registerTIADataStore
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/25/04
//-----------------------------------------------------------------------------
bool N_TIA_MPDEInterface::registerTIADataStore(N_TIA_DataStore * ds_tmp)

{
  dsPtr_ = ds_tmp;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::registerTIAControl
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/25/04
//-----------------------------------------------------------------------------
bool N_TIA_MPDEInterface::registerTIAControl(N_ANP_AnalysisManager * anaManager_tmp)

{
  anaManagerPtr_ = anaManager_tmp;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::registerTIAStepErrorControl
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 1414
// Creation Date : 10/04/05
//-----------------------------------------------------------------------------
  bool N_TIA_MPDEInterface::registerTIAStepErrorControl(N_TIA_StepErrorControl * tiaSec_tmp)

{
  tiaSecPtr_ = tiaSec_tmp;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::setInitialCondition
// Purpose       : Method to specify transient initial condition
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/25/04
//-----------------------------------------------------------------------------
bool N_TIA_MPDEInterface::setInitialCondition(N_LAS_Vector * initialConditionPtr)
{
  *(dsPtr_->nextSolutionPtr) = *(initialConditionPtr);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::setStateInitialCondition
// Purpose       : Method to specify transient state initial condition
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
bool N_TIA_MPDEInterface::setStateInitialCondition(N_LAS_Vector * stateInitialConditionPtr)
{
  *(dsPtr_->nextStatePtr) = *(stateInitialConditionPtr);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::setStoreInitialCondition
// Purpose       : Method to specify transient store initial condition
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool N_TIA_MPDEInterface::setStoreInitialCondition(N_LAS_Vector * storeInitialConditionPtr)
{
  *(dsPtr_->nextStorePtr) = *(storeInitialConditionPtr);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::setInitialCondition
// Purpose       : Method to specify transient initial condition
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/25/04
//-----------------------------------------------------------------------------
bool N_TIA_MPDEInterface::setQVectorInitialCondition(N_LAS_Vector * qVectorInitialConditionPtr)
{
  *(dsPtr_->daeQVectorPtr) = *(qVectorInitialConditionPtr);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::runDCOP
// Purpose       : Execute the top level DC sweep control loop.
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/25/04
//-----------------------------------------------------------------------------
bool N_TIA_MPDEInterface::runDCOP()
{
  if( !tiaParams_.NOOP )   
  {
    // as long as the user didn't request "noop"
    // try and do the operating point calculation
    anaManagerPtr_->analysis = ANP_MODE_DC_SWEEP;
    return anaManagerPtr_->run();
  }
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::runTransient
// Purpose       : Execute the top level transient control loop without DCOP.
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/25/04
//-----------------------------------------------------------------------------
bool N_TIA_MPDEInterface::runTransient()
{
  anaManagerPtr_->analysis = ANP_MODE_TRANSIENT;
  tiaParams_.NOOP = true;
  return anaManagerPtr_->run();
}


//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::runTransient
// Purpose       : Execute the top level transient control loop with DCOP.
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/25/04
//-----------------------------------------------------------------------------
bool N_TIA_MPDEInterface::runTransientWithDCOP()
{
  anaManagerPtr_->analysis = ANP_MODE_TRANSIENT;
  return anaManagerPtr_->run();
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::getTimeIntInfo
// Purpose       : Get details from time integrator on last step.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, 1437
// Creation Date : 10/22/07
//-----------------------------------------------------------------------------
void N_TIA_MPDEInterface::getTimeIntInfo(N_TIA_TimeIntInfo & tiInfo)
{
  anaManagerPtr_->getTimeIntInfo(tiInfo);
}


//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::runStep
// Purpose       : Take one, prescribed, time-step
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, 1437
// Creation Date : 10/22/07
//-----------------------------------------------------------------------------
bool N_TIA_MPDEInterface::runStep(const N_TIA_TimeIntInfo & tiInfo, N_TIA_TwoLevelError & tlError)
{
  anaManagerPtr_->startTimeStep(tiInfo);
  return anaManagerPtr_->runStep(tiInfo, tlError);
}


//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::stepSuccess
// Purpose       : Process a successful step
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, 1437
// Creation Date : 10/22/07
//-----------------------------------------------------------------------------
void N_TIA_MPDEInterface::stepSuccess(int analysisType)
{
  anaManagerPtr_->stepSuccess( analysisType );
}


//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::getFinalSolution
// Purpose       : Method to return output solution (either from DCOP or transient)
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/25/04
//-----------------------------------------------------------------------------
N_LAS_Vector* N_TIA_MPDEInterface::getFinalSolution()
{
  return dsPtr_->currSolutionPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::getStateFinalSolution
// Purpose       : Method to return output state vector (either from DCOP or transient)
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
N_LAS_Vector* N_TIA_MPDEInterface::getStateFinalSolution()
{
  return dsPtr_->currStatePtr;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::getStoreFinalSolution
// Purpose       : Method to return output store vector (either from DCOP or transient)
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
N_LAS_Vector* N_TIA_MPDEInterface::getStoreFinalSolution()
{
  return dsPtr_->currStorePtr;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::getQVectorFinalSolution
// Purpose       : Method to return output Q vector (either from DCOP or transient)
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
N_LAS_Vector* N_TIA_MPDEInterface::getQVectorFinalSolution()
{
  return dsPtr_->daeQVectorPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::getQVectorHistoryFinalSolution
// Purpose       : Method to return output Q vector (either from DCOP or transient)
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
N_LAS_Vector* N_TIA_MPDEInterface::getQVectorHistoryFinalSolution()
{
  return dsPtr_->qHistory[0];
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::reInitialize
// Purpose       : Method to re-initialize time integrator in between separate
//               : integrations with the same size vectors.
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/25/04
//-----------------------------------------------------------------------------
bool N_TIA_MPDEInterface::reInitialize()
{
  tiaParams_.resume = false;
  anaManagerPtr_->setBeginningIntegrationFlag(true);
  // 10/04/05 tscofffe:  We need to reinitialize the variables in the
  // StepErrorControl object.
#ifdef Xyce_DEBUG_TIME
  cout << *(tiaSecPtr_) << endl;
#endif
  tiaSecPtr_->resetAll();
#ifdef Xyce_DEBUG_TIME
  cout << *(tiaSecPtr_) << endl;
#endif

  // 10/04/05 tscoffe:  BackwardDifferentiation15 needs to be re-initialized here.
  anaManagerPtr_->wimPtr->initialize();
  return true;
}

