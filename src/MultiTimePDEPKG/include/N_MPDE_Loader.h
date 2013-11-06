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
// Filename       : $RCSfile: N_MPDE_Loader.h,v $
//
// Purpose        : MPDE Specific Loader
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/11/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.32.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:46 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_MPDE_Loader_H
#define Xyce_MPDE_Loader_H

// ---------- Standard Includes ----------

#include <vector>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>
#include <N_LOA_Loader.h>

#include <N_MPDE_State.h>
#include <N_MPDE_WarpedPhaseCondition.h>
#include <N_MPDE_DeviceInterface.h>

// ---------- Forward declarations --------

class N_LAS_Vector;
class N_LAS_BlockVector;
class N_LAS_Matrix;
class N_LAS_BlockMatrix;

class N_MPDE_Discretization;
class N_MPDE_Manager;

//-----------------------------------------------------------------------------
// Class         : N_MPDE_Loader
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
class N_MPDE_Loader : public N_LOA_Loader
{
public:

  // Default constructor
  N_MPDE_Loader( N_MPDE_State & state,
                 const RefCountPtr<const N_MPDE_Discretization> disc,
                 const RefCountPtr<N_MPDE_Manager> mgr,
                 const Teuchos::RefCountPtr<N_MPDE_WarpedPhaseCondition>& warpMPDEPhasePtr)
  : state_(state),
    fastTimeDiscRCPtr_(disc),
    periodicTimesOffset_(0),
    period_(0),
    bOmegadQdt2Ptr_(0),
    mpdeMgrRCPtr_(mgr),
    warpMPDEPhasePtr_(warpMPDEPhasePtr),
    warpMPDE_(warpMPDEPhasePtr != Teuchos::null)
  {}

  // Destructor
  ~N_MPDE_Loader();

  // Method which is called to load the new-DAE contributions 
  bool loadDAEMatrices( N_LAS_Vector * X,
                        N_LAS_Vector * S,
                        N_LAS_Vector * dSdt,
                        N_LAS_Vector * Store,
                        N_LAS_Matrix * dQdx,
                        N_LAS_Matrix * dFdx );

  // Method which is called to load the new-DAE vectors
  bool loadDAEVectors( N_LAS_Vector * X,
                       N_LAS_Vector * currX,
                       N_LAS_Vector * lastX,
                       N_LAS_Vector * S,
                       N_LAS_Vector * currS,
                       N_LAS_Vector * lastS,
                       N_LAS_Vector * dSdt,
                       N_LAS_Vector * Store,
                       N_LAS_Vector * currStore,
                       N_LAS_Vector * lastStore,
                       N_LAS_Vector * storeLeadCurrQComp,
                       N_LAS_Vector * Q,
                       N_LAS_Vector * F,
                       N_LAS_Vector * dFdxdVp,
                       N_LAS_Vector * dQdxdVp );

  bool updateState(    N_LAS_Vector * nextSolVectorPtr,
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
  void setFastTimes( const vector<double> & times );
  
  void setPeriodFlags( const vector<bool> & periodicFlags );
  
  // Construct a periodic version of times_
  void constructPeriodicTimes();

  // Registration method for the device packaage
  void registerAppLoader( RefCountPtr<N_LOA_Loader> appLoaderPtr );
  
  // Register the interface to the devices
  void registerMPDEDeviceInterface( RefCountPtr< N_MPDE_DeviceInterface > mpdeDevIntPtr );

  // Registration method for the next Application Vector
  void registerAppNextVec( RefCountPtr<N_LAS_Vector> appNextVecPtr );

  // Registration method for the curr Application Vector
  void registerAppCurrVec( RefCountPtr<N_LAS_Vector> appCurrVecPtr );

  // Registration method for the last Application Vector
  void registerAppLastVec( RefCountPtr<N_LAS_Vector> appLastVecPtr );

  // Registration method for the next Application State Vector
  void registerAppNextStaVec( RefCountPtr<N_LAS_Vector> appNextStaVecPtr );

  // Registration method for the curr Application State Vector
  void registerAppCurrStaVec( RefCountPtr<N_LAS_Vector> appCurrStaVecPtr );

  // Registration method for the last Application State Vector
  void registerAppLastStaVec( RefCountPtr<N_LAS_Vector> appLastStaVecPtr );

  // Registration method for the next Application Store Vector
  void registerAppNextStoVec( RefCountPtr<N_LAS_Vector> appNextStoVecPtr );

  // Registration method for the curr Application Store Vector
  void registerAppCurrStoVec( RefCountPtr<N_LAS_Vector> appCurrStoVecPtr );

  // Registration method for the curr Application Store Vector
  void registerAppLastStoVec( RefCountPtr<N_LAS_Vector> appLastStoVecPtr );
  
  // Registration method for the Application Store Vector Lead Current Q Component 
  void registerAppStoLeadCurrQCompVec( RefCountPtr<N_LAS_Vector> appStoLeadCurrQCompVecPtr );
  
  // Registration method for the Application dQdx Matrix
  void registerAppdQdx( RefCountPtr<N_LAS_Matrix> appdQdxPtr );

  // Registration method for the Application dFdx Matrix
  void registerAppdFdx( RefCountPtr<N_LAS_Matrix> appdFdxPtr );

  // Registration method for the MPDE size dQdx Block Matrix
  void registerMPDEdQdx( RefCountPtr<N_LAS_BlockMatrix> bmdQdxPtr );

  // Registration method for the MPDE size dFdx Block Matrix
  void registerMPDEdFdx( RefCountPtr<N_LAS_BlockMatrix> bmdFdxPtr );

  // Registration method for a temp vector in constructin dq/dt2
  void registerOmegadQdt2( RefCountPtr<N_LAS_BlockVector> omegadQdt2Ptr);

private :

  //MPDE State
  N_MPDE_State & state_;

  // discretization
  RefCountPtr<const N_MPDE_Discretization> fastTimeDiscRCPtr_;

  //Fast Time Scale Points
  vector<double> times_;
  int periodicTimesOffset_;
  vector<double> periodicTimes_;
  double period_;
  
  //a vector of bools to indicate if a solution variable does not
  //appear to be periodic.  Non-periodic signals will be handled
  //differently by the MPDE_Loader class
  vector<bool> nonPeriodic_;
  
  // Base Application loader
  RefCountPtr<N_LOA_Loader> appLoaderPtr_;
  RefCountPtr < N_MPDE_DeviceInterface > mpdeDevInterfacePtr_;

  // Application Linear Objects
  RefCountPtr<N_LAS_Vector> appNextVecPtr_;
  RefCountPtr<N_LAS_Vector> appCurrVecPtr_;
  RefCountPtr<N_LAS_Vector> appLastVecPtr_;

  RefCountPtr<N_LAS_Vector> appNextStaVecPtr_;
  RefCountPtr<N_LAS_Vector> appCurrStaVecPtr_;
  RefCountPtr<N_LAS_Vector> appLastStaVecPtr_;
  RefCountPtr<N_LAS_Matrix> appdQdxPtr_;
  RefCountPtr<N_LAS_Matrix> appdFdxPtr_;
  RefCountPtr<N_LAS_Vector> appNextStoVecPtr_;
  RefCountPtr<N_LAS_Vector> appCurrStoVecPtr_;
  RefCountPtr<N_LAS_Vector> appLastStoVecPtr_;
  Teuchos::RCP<N_LAS_Vector> appStoLeadCurrQCompVecPtr_;

  // Tmp storage block matrices 
  RefCountPtr<N_LAS_BlockMatrix> bmdQdxPtr_;
  RefCountPtr<N_LAS_BlockMatrix> bmdFdxPtr_;

  // mpde manager:
  RefCountPtr<N_MPDE_Manager> mpdeMgrRCPtr_;

  // MPDE/WaMPDE data
  RefCountPtr<N_LAS_BlockVector> bOmegadQdt2Ptr_;
  
  bool warpMPDE_;
  Teuchos::RefCountPtr<N_MPDE_WarpedPhaseCondition> warpMPDEPhasePtr_;
  
  //bool N_MPDE_Loader::findNonPeriodicSignals_(const N_LAS_BlockVector & solutionBlockVector );

};


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppLoader
// Purpose       : Registration method for the device packaage
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppLoader( RefCountPtr<N_LOA_Loader> appLoaderPtr )
{ 
  appLoaderPtr_ = appLoaderPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerMPDEDeviceInterface
// Purpose       : Registration method for the device packaage
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerMPDEDeviceInterface( RefCountPtr< N_MPDE_DeviceInterface > mpdeDevIntPtr )
{ 
  mpdeDevInterfacePtr_ = mpdeDevIntPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppNextVec
// Purpose       : Registration method for the Application Vector
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader:: registerAppNextVec( RefCountPtr<N_LAS_Vector> appNextVecPtr )
{ 
  appNextVecPtr_ = appNextVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppCurrVec
// Purpose       : Registration method for the Application Vector
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 1437
// Creation Date : 
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader:: registerAppCurrVec( RefCountPtr<N_LAS_Vector> appCurrVecPtr )
{ 
  appCurrVecPtr_ = appCurrVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppLastVec
// Purpose       : Registration method for the Application Vector
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 1437
// Creation Date :
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader:: registerAppLastVec( RefCountPtr<N_LAS_Vector> appLastVecPtr )
{ 
  appLastVecPtr_ = appLastVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppNextStaVec
// Purpose       : Registration method for the Application State Vector
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppNextStaVec( RefCountPtr<N_LAS_Vector> appNextStaVecPtr )
{ 
  appNextStaVecPtr_ = appNextStaVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppCurrStaVec
// Purpose       : Registration method for the Application State Vector
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 01/25/07
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppCurrStaVec( RefCountPtr<N_LAS_Vector> appCurrStaVecPtr )
{ 
  appCurrStaVecPtr_ = appCurrStaVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppLastStaVec
// Purpose       : Registration method for the Application State Vector
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 01/25/07
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppLastStaVec( RefCountPtr<N_LAS_Vector> appLastStaVecPtr )
{ 
  appLastStaVecPtr_ = appLastStaVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppNextStoVec
// Purpose       : Registration method for the Application Store Vector
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date :
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppNextStoVec( RefCountPtr<N_LAS_Vector> appNextStoVecPtr )
{ 
  appNextStoVecPtr_ = appNextStoVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppCurrStoVec
// Purpose       : Registration method for the Application Store Vector
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppCurrStoVec( RefCountPtr<N_LAS_Vector> appCurrStoVecPtr )
{ 
  appCurrStoVecPtr_ = appCurrStoVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppLastStoVec
// Purpose       : Registration method for the Application Store Vector
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppLastStoVec( RefCountPtr<N_LAS_Vector> appLastStoVecPtr )
{ 
  appLastStoVecPtr_ = appLastStoVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppStoLeadCurrQCompVec
// Purpose       : Registration method for the Application Store Vector
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek
// Creation Date : 
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppStoLeadCurrQCompVec( RefCountPtr<N_LAS_Vector> appStoLeadCurrQCompVecPtr )
{ 
  appStoLeadCurrQCompVecPtr_ = appStoLeadCurrQCompVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerOmegadQdt2
// Purpose       : Registration method for a temp vector in constructin dq/dt2
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414, Computational Sciences
// Creation Date : 08/10/05
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerOmegadQdt2( RefCountPtr<N_LAS_BlockVector> omegadQdt2Ptr)
{ 
  bOmegadQdt2Ptr_ = omegadQdt2Ptr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppdQdx
// Purpose       : Registration method for the Application dQdx Matrix
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppdQdx( RefCountPtr<N_LAS_Matrix> appdQdxPtr )
{ 
  appdQdxPtr_ = appdQdxPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppdFdx
// Purpose       :  Registration method for the Application dFdx Matrix
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppdFdx( RefCountPtr<N_LAS_Matrix> appdFdxPtr )
{ 
  appdFdxPtr_ = appdFdxPtr; 
}

// Registration method for the MPDE size dQdx Block Matrix
inline void N_MPDE_Loader::registerMPDEdQdx( RefCountPtr<N_LAS_BlockMatrix> bmdQdxPtr )
{ 
  bmdQdxPtr_ = bmdQdxPtr; 
}

// Registration method for the MPDE size dFdx Block Matrix
inline void N_MPDE_Loader::registerMPDEdFdx( RefCountPtr<N_LAS_BlockMatrix> bmdFdxPtr )
{ 
  bmdFdxPtr_ = bmdFdxPtr; 
}

#endif

