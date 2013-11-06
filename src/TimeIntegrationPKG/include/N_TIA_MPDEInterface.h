
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
// Filename      : $RCSfile: N_TIA_MPDEInterface.h,v $
//
// Purpose       : This file defines the time integration interface for MPDE
//
// Special Notes :
//
// Creator       : Todd S. Coffey, 9214
//
// Creation Date : 03/24/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.19.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_TIME_MPDEINTERFACE_H
#define Xyce_TIME_MPDEINTERFACE_H

// ---------- Standard Includes ----------

#include <list>
#include <vector>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>

#include <N_ANP_AnalysisManager.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_TimeIntegrationMethods.h>

#include <N_LAS_Vector.h>

#include <N_IO_PkgOptionsMgr.h>

// ---------- Forward Declarations ----------

class N_TIA_TIAParams;

class N_ANP_AnalysisManager;

class N_TIA_StepErrorControl;
class N_TIA_DataStore;
class N_TIA_TwoLevelError;
class N_TIA_TimeIntInfo;

//-----------------------------------------------------------------------------
// Class         : N_TIA_MPDEInterface
// Purpose       : This is the interface for the time integrator that is
//               : provided to the MPDE loader.
// Special Notes :
// Creator       : Todd S. Coffey, 9214
// Creation Date : 03/24/04
//-----------------------------------------------------------------------------
class N_TIA_MPDEInterface
{
public:

  // Default donstructor
  N_TIA_MPDEInterface(N_TIA_TIAParams & tp);

  // Destructor
  ~N_TIA_MPDEInterface();

  // Execution functions:

  // Method to register the TIA parameters object.
  bool registerTIAParams(const N_TIA_TIAParams & tiaParams_tmp);

  // Method to register the TIA Data Store object.
  bool registerTIADataStore(N_TIA_DataStore * ds_tmp);

  // Method to register the TIA Control Alogrithm object.
  bool registerTIAControl(N_ANP_AnalysisManager * tiaControl_tmp);
  //
  // Method to register the TIA StepErrorControl object.
  bool registerTIAStepErrorControl(N_TIA_StepErrorControl * tiaSec_tmp);

  // Method to specify transient initial condition
  bool setInitialCondition(N_LAS_Vector * initialConditionPtr);

  // Method to specify transient state initial condition
  bool setStateInitialCondition(N_LAS_Vector *stateInitialConditionPtr);

  // Method to specify transient store initial condition
  bool setStoreInitialCondition(N_LAS_Vector *storeInitialConditionPtr);

  // Method to specify transient Q initial condition
  bool setQVectorInitialCondition(N_LAS_Vector *qVectorInitialConditionPtr);
  
  // Method to run DCOP
  bool runDCOP();

  // Method to run transient without DCOP
  bool runTransient();
  
  // Method to run transient with DCOP
  bool runTransientWithDCOP();

  // Methods to take single time steps
  void getTimeIntInfo(N_TIA_TimeIntInfo & tiInfo);
  bool runStep(const N_TIA_TimeIntInfo & tiInfo, N_TIA_TwoLevelError & tlError);
  void stepSuccess(int analysis_type);
  
  // Method to return output solution (either from DCOP or transient)
  N_LAS_Vector* getFinalSolution();

  // Method to return output state vector (either from DCOP or transient)
  N_LAS_Vector* getStateFinalSolution();

  // Method to return output store vector (either from DCOP or transient)
  N_LAS_Vector* getStoreFinalSolution();

  // Method to return output DAE Q vector
  N_LAS_Vector* getQVectorFinalSolution();
  
  // Method to return output DAE Q vector history
  N_LAS_Vector* getQVectorHistoryFinalSolution();
  
  // Method to re-initialize time integrator
  bool reInitialize();

  // Accessor to the N_TIA_TIAParams object.
  N_TIA_TIAParams& getTIAParams();

  bool setTIAParams(N_TIA_TIAParams& tiaParams);

protected:

private :

  // Pointer to DataStore
  N_TIA_DataStore * dsPtr_;  

  // Pointer to TIAParams
  N_TIA_TIAParams & tiaParams_;

  // Pointer to the Analysis Manager
  N_ANP_AnalysisManager * anaManagerPtr_;
  
  // Pointer to the TIA StepErrorControl object
  N_TIA_StepErrorControl * tiaSecPtr_;

};

//-----------------------------------------------------------------------------
// Function      : N_TIA_MPDEInterface::getTIAParams
// Purpose       : Accessor to the N_TIA_TIAParams object.
// Special Notes : 
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 3/25/04
//-----------------------------------------------------------------------------
inline N_TIA_TIAParams& N_TIA_MPDEInterface::getTIAParams()
{
  return (tiaParams_);
}

#endif // Xyce_TIME_MPDEINTERFACE_H

