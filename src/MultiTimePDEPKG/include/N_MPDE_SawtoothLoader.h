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
// Filename       : $RCSfile: N_MPDE_SawtoothLoader.h,v $
//
// Purpose        : MPDE Sawtooth IC Specific Loader
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/23/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.18.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:46 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_MPDE_SawtoothLoader_H
#define Xyce_MPDE_SawtoothLoader_H

// ---------- Standard Includes ----------

#include <vector>
#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>
#include <N_LOA_Loader.h>

#include <N_MPDE_State.h>

// ---------- Forward declarations --------

class N_LAS_Vector;
class N_LAS_Matrix;

class N_ANP_AnalysisInterface;

//-----------------------------------------------------------------------------
// Class         : N_MPDE_SawtoothLoader
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/23/04
//-----------------------------------------------------------------------------
class N_MPDE_SawtoothLoader : public N_LOA_Loader
{
public:

  // Default constructor
  N_MPDE_SawtoothLoader( N_MPDE_State & state )
  : state_(state),
    appLoader_(0),
    anaInt_(0),
    timeShift_(0.0)
  {}

  // Destructor
  ~N_MPDE_SawtoothLoader() {}

  // Method which is called to load the new-DAE contributions 
  bool loadDAEMatrices( N_LAS_Vector * X,
                        N_LAS_Vector * S,
                        N_LAS_Vector * dSdt,
                        N_LAS_Vector * Store,
                        N_LAS_Matrix * dQdx,
                        N_LAS_Matrix * dFdx);

  // Method which is called to load the new-DAE vectors
  bool loadDAEVectors( N_LAS_Vector * X,
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
                       N_LAS_Vector * dQdxdVp );

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

  // Virtual method which initializes the nonlinear problem.
  virtual bool initializeProblem( N_LAS_Vector * nextSolVectorPtr,
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
                          N_LAS_Vector * dQdxdVpVectorPtr)  { return false; }

  // Assign times for fast time scale
  void setTimeShift( double timeShift )
  { timeShift_ = timeShift; }

  // Registration method for the system loader 
  void registerAppLoader( RefCountPtr<N_LOA_Loader> appLoader )
  { appLoader_ = appLoader; }

  // Registration method for the TIA
  void registerTIA( RefCountPtr<N_ANP_AnalysisInterface> anaInt )
  { anaInt_ = anaInt; }

private :

  //MPDE State
  N_MPDE_State & state_;

  //Time Shift for Fast Src
  double timeShift_;

  // Base Application loader
  RefCountPtr<N_LOA_Loader> appLoader_;

  // TIA Manager
  RefCountPtr<N_ANP_AnalysisInterface> anaInt_;

};

#endif

