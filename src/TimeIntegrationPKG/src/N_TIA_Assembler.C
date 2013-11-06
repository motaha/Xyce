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
// Filename      : $RCSfile: N_TIA_Assembler.C,v $
// Purpose       : 
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 1/1/07
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.24.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:49 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_TIA_Assembler.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_TIA_DataStore.h>

#include <N_LOA_Loader.h>

#include <N_UTL_Timer.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

//-----------------------------------------------------------------------------
// Function      : N_TIA_Assembler::N_TIA_Assembler::
// Purpose       : constructor
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 01/30/07
//-----------------------------------------------------------------------------
N_TIA_Assembler::N_TIA_Assembler ( N_TIA_DataStore & ds,
                                   N_LOA_Loader & loader,
                                   N_TIA_WorkingIntegrationMethod & wim,
                                   N_PDS_Manager & pds,
                                   bool daeStateDerivFlag
                                   )
: ds_(ds),
  loader_(loader),
  wim_(wim),
  pdsMgr(pds),
  residualTimerPtr_(0),
  jacobianTimerPtr_(0),
  daeStateDerivFlag_(daeStateDerivFlag)
{
  residualTimerPtr_ = new N_UTL_Timer(*(pdsMgr.getPDSComm()));
  jacobianTimerPtr_ = new N_UTL_Timer(*(pdsMgr.getPDSComm()));
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_Assembler::~N_TIA_Assembler::
// Purpose       : destructor
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 01/30/07
//-----------------------------------------------------------------------------

N_TIA_Assembler::~N_TIA_Assembler()
{
  if (residualTimerPtr_ != 0)  delete residualTimerPtr_;
  if (jacobianTimerPtr_ != 0)  delete jacobianTimerPtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DAE_Assembler::N_TIA_DAE_Assembler::
// Purpose       : constructor
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 01/30/07
//-----------------------------------------------------------------------------
N_TIA_DAE_Assembler::N_TIA_DAE_Assembler 
( N_TIA_DataStore & ds,
  N_LOA_Loader & loader,
  N_TIA_WorkingIntegrationMethod & wim,
  N_PDS_Manager & pds,
  bool daeStateDerivFlag
  )
 : N_TIA_Assembler(ds,loader,wim,pds,daeStateDerivFlag)
{

}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DAE_Assembler::~N_TIA_DAE_Assembler::
// Purpose       : destructor
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 01/30/07
//-----------------------------------------------------------------------------

N_TIA_DAE_Assembler::~N_TIA_DAE_Assembler()
{
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DAE_Assembler::loadRHS
//
// Purpose       : This function manages the various function calls necessary
//                 to assemble the new-DAE form of the residual (RHS).
//
// Special Notes : All the contributions to the RHS come (in the new DAE
//                 form) from the device package, as Q, F, and B.  The 
//                 RHS needs dQdt + F - B.  As dQdt is determined by the
//                 time integration package, the final summation should be
//                 managed from here.
//
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 03/04/04
//-----------------------------------------------------------------------------
bool N_TIA_DAE_Assembler::loadRHS ()
{
  bool bsuccess = true;
  bool tmpBool = true;

  // Start the timer...
  residualTimerPtr_->resetStartTime();

  ds_.daeQVectorPtr->putScalar(0.0);
  ds_.daeFVectorPtr->putScalar(0.0);

  ds_.dFdxdVpVectorPtr->putScalar(0.0);
  ds_.dQdxdVpVectorPtr->putScalar(0.0);
  ds_.dQdxMatrixPtr->put(0.0);
  ds_.dFdxMatrixPtr->put(0.0);

  // Update the state. Note - underneath this call, most of the calculations
  // pertaining to the currents, conductances, etc. will happen.
  tmpBool = loader_.updateState 
              ((ds_.nextSolutionPtr), 
               (ds_.currSolutionPtr), 
               (ds_.lastSolutionPtr), 
               (ds_.nextStatePtr),
               (ds_.currStatePtr),
               (ds_.lastStatePtr),
               (ds_.nextStorePtr),
               (ds_.currStorePtr),
               (ds_.lastStorePtr)
               );

  bsuccess = bsuccess && tmpBool;

  if (daeStateDerivFlag_)
  {
    wim_.updateStateDeriv ();
  }

  // first load the 2 components: Q  and F  
  // (there is no longer a separate B vector for indep sources)
  tmpBool = loader_.loadDAEVectors   
              ((ds_.nextSolutionPtr), 
               (ds_.currSolutionPtr), 
               (ds_.lastSolutionPtr), 
               (ds_.nextStatePtr),
               (ds_.currStatePtr),
               (ds_.lastStatePtr),
               (ds_.nextStateDerivPtr),
               (ds_.nextStorePtr),
               (ds_.currStorePtr),
               (ds_.lastStorePtr),
               (ds_.nextStoreLeadCurrQCompPtr),
               (ds_.daeQVectorPtr),
               (ds_.daeFVectorPtr),
               (ds_.dFdxdVpVectorPtr),
               (ds_.dQdxdVpVectorPtr) );
  bsuccess = bsuccess && tmpBool;
  wim_.updateLeadCurrent();
  // Now determine dQdt:
  // now sum them all together, to create the total.
  // f(x) is given by:
  //
  //    f(x) = dQ/dt + F(x) = 0
  //
  // Note, the nonlinear solver is expecting the RHS vector to
  // contain -f(x).  Or, possibly -f(x) + J*dx, if voltage 
  // limiting is on.
  wim_.obtainResidual();

  // Update the total load time
  residualTime_ = residualTimerPtr_->elapsedTime();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DAE_Assembler::loadJacobian
//
// Purpose       : This function manages the various function calls necessary
//                 to assemble the new-DAE form of the Jacobian.
//
// Special Notes : All the contributions to the Jacobian 
//                 come (in the new DAE form)
//                 from the device package, as dQdx and dFdx. 
//
//                 The Jacobian is: 
//
//                 J = df/dx = d(dQdt)/dx + dF/dx 
//
//                 As dQdt is determined by the time integration package, 
//                 the final summation should be managed from here.
//
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 03/04/04
//-----------------------------------------------------------------------------
bool N_TIA_DAE_Assembler::loadJacobian ()
{
  bool bsuccess = true;
  bool tmpBool = true;

  // Start the timer...
  jacobianTimerPtr_->resetStartTime();

#if 0
  // This put(0) isn't needed anymore, as it is handled down in the 
  // obtainJacobian method.
  ds_.JMatrixPtr->put(0.0);  
#endif

#if 0  
 // don't put these here!  In the MPDE case they are loaded 
 // during the vector load.
  ds_.dQdxMatrixPtr->put(0.0);
  ds_.dFdxMatrixPtr->put(0.0);
#endif

  // first load the 2 components: dQdx and dFdx
  tmpBool = loader_.loadDAEMatrices 
              ((ds_.nextSolutionPtr),
               (ds_.nextStatePtr),
               (ds_.nextStateDerivPtr),
               (ds_.nextStorePtr),
               (ds_.dQdxMatrixPtr),
               (ds_.dFdxMatrixPtr));
  bsuccess = bsuccess && tmpBool;

  // Now determine the d(dQdt)/dx stuff:

  // now sum them all together, to create the total:
  // J = alpha/dt * dQ/dx + dF/dx
  wim_.obtainJacobian();

  // Update the total load time
  jacobianTime_ = jacobianTimerPtr_->elapsedTime();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DAE_Assembler::applyJacobian
//
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
bool N_TIA_DAE_Assembler::applyJacobian (const N_LAS_Vector& input, N_LAS_Vector& result)
{
  bool bsuccess = true;
  bool tmpBool = true;

  // Start the timer...
  jacobianTimerPtr_->resetStartTime();

  // first load the 2 components: dQdx and dFdx
  tmpBool = loader_.applyDAEMatrices 
              ((ds_.nextSolutionPtr),
               (ds_.nextStatePtr),
               (ds_.nextStateDerivPtr),
               (ds_.nextStorePtr),
               input,
               (ds_.dQdxVecVectorPtr),
               (ds_.dFdxVecVectorPtr));
  bsuccess = bsuccess && tmpBool;

  // Now determine the d(dQdt)/dx stuff:

  // now sum them all together, to create the total:
  // J = alpha/dt * dQ/dx + dF/dx
  wim_.applyJacobian(input, result);

  // Update the total load time
  jacobianTime_ = jacobianTimerPtr_->elapsedTime();

  return bsuccess;
}

