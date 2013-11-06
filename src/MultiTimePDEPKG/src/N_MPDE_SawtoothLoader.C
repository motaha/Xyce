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
// Filename      : $RCSfile: N_MPDE_SawtoothLoader.C,v $
// Purpose       : 
// Special Notes :
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/23/04
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.21.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:47 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <iostream>

// ----------   Xyce Includes   ----------

#include <N_MPDE_SawtoothLoader.h>

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_ANP_AnalysisInterface.h>

//-----------------------------------------------------------------------------
// Function      : N_MPDE_SawtoothLoader::loadDAEMatrices
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/23/04
//-----------------------------------------------------------------------------
bool N_MPDE_SawtoothLoader::loadDAEMatrices( N_LAS_Vector * X,
                                             N_LAS_Vector * S,
                                             N_LAS_Vector * dSdt,
                                             N_LAS_Vector * Store,
                                             N_LAS_Matrix * dQdx,
                                             N_LAS_Matrix * dFdx)
{
  //Zero out matrices
  dQdx->put(0.0);
  dFdx->put(0.0);

  double fastTime = anaInt_->getTime() + timeShift_;
  //Set Time for fast time scale somewhere
  state_.fastTime = fastTime;

  return appLoader_->loadDAEMatrices( X, S, dSdt, Store, dQdx, dFdx);
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_SawtoothLoader::loadDAEVectors
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_MPDE_SawtoothLoader::loadDAEVectors( N_LAS_Vector * X,
                                            N_LAS_Vector * currX,
                                            N_LAS_Vector * lastX,
                                            N_LAS_Vector * nextS,
                                            N_LAS_Vector * currS,
                                            N_LAS_Vector * lastS,
                                            N_LAS_Vector * dSdt,
                                            N_LAS_Vector * nextStore,
                                            N_LAS_Vector * currStore,
                                            N_LAS_Vector * lastStore,
                                            N_LAS_Vector * storeLeadCurrQComp,
                                            N_LAS_Vector * Q,
                                            N_LAS_Vector * F,
                                            N_LAS_Vector * dFdxdVp,
                                            N_LAS_Vector * dQdxdVp )
{
  double fastTime = anaInt_->getTime() + timeShift_;
  //Set Time for fast time scale somewhere
  state_.fastTime = fastTime;

  appLoader_->updateSources();

  return appLoader_->loadDAEVectors
    ( X, currX, lastX, nextS, currS, lastS, dSdt, 
      nextStore, currStore, lastStore, storeLeadCurrQComp, Q, F, dFdxdVp, dFdxdVp );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_SawtoothLoader::updateState
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool N_MPDE_SawtoothLoader::updateState (N_LAS_Vector * nextSolVectorPtr,
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
  return appLoader_->updateState
    ( nextSolVectorPtr, currSolVectorPtr, lastSolVectorPtr,
      nextStaVectorPtr, currStaVectorPtr, lastStaVectorPtr,
      nextStoVectorPtr, currStoVectorPtr, lastStoVectorPtr);
}

