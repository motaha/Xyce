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
// Filename       : $RCSfile: N_LOA_HBLoader.h,v $
//
// Purpose        : HB Specific Loader
//
// Special Notes  :
//
// Creator        : Todd Coffey, Ting Mei
//
// Creation Date  : 07/28/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:46 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_LOA_HBLoader_H
#define Xyce_LOA_HBLoader_H

// ---------- Standard Includes ----------

#include <vector>

#include <Teuchos_RCP.hpp>

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>
#include <N_LOA_Loader.h>

#include <N_MPDE_State.h>
#include <N_DEV_DeviceInterface.h>

// ---------- Forward declarations --------

class N_LAS_Vector;
class N_LAS_BlockVector;
class N_LAS_Matrix;
class N_LAS_BlockMatrix;

class N_MPDE_Discretization;

class N_LAS_HBBuilder;

//-----------------------------------------------------------------------------
// Class         : N_LOA_HBLoader
// Purpose       : HB specific CktLoader interface
// Special Notes :
// Creator       : Todd Coffey, Ting Mei, Rich Schiek
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
class N_LOA_HBLoader : public N_LOA_Loader
{
public:

  // Default constructor
    N_LOA_HBLoader( N_MPDE_State & state,
                 const Teuchos::RCP<const N_MPDE_Discretization> discRCPtr )
  : state_(state),
    fastTimeDiscRCPtr_(discRCPtr),
    periodicTimesOffset_(0),
    period_(0),
    matrixFreeFlag_(false)
  {}

  // Destructor
  ~N_LOA_HBLoader();

  // Method which is called to load the new-DAE contributions 
  bool loadDAEMatrices( N_LAS_Vector * X,
                        N_LAS_Vector * S,
                        N_LAS_Vector * dSdt,
                        N_LAS_Vector * Store,
                        N_LAS_Matrix * dQdx,
                        N_LAS_Matrix * dFdx );

  // Method which is called to load the time-dependent HB Matrices
  bool loadTimeDepDAEMatrices( N_LAS_Vector * X,
                        N_LAS_Vector * S,
                        N_LAS_Vector * dSdt,
                        N_LAS_Vector * Store,
                        N_LAS_Matrix * dQdx,
                        N_LAS_Matrix * dFdx );

  // Method for matrix-free 
  bool applyDAEMatrices( N_LAS_Vector * X,
                         N_LAS_Vector * S,
                         N_LAS_Vector * dSdt,
                         N_LAS_Vector * Store,
                         const N_LAS_Vector & V,
                         N_LAS_Vector * dQdxV,
                         N_LAS_Vector * dFdxV );

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

  // New functions for LOA_HBLoader
  // Assumption:  xt is in the block format:
  // xt.block(i) = { x_n(t_i) }_{n=0..N} N = number of solution components, t_i = time point i
  // P takes xt's block format and converts it to: 
  // (P*xt).block(n) = { x_n(t_i) }_{i=0..T} T = number of time points in block vector.
 
  void setFastTimes( const vector<double> & times );
  void setPeriodFlags( const vector<bool> & periodicFlags );
  void constructPeriodicTimes();

  // Is the problem being applied as a matrix free operator.
  bool matrixFreeFlag() const { return matrixFreeFlag_; }
  void setMatrixFreeFlag(bool matrixFreeFlag) { matrixFreeFlag_ = matrixFreeFlag; }
 
  // xf = D*P*xt, xf has the same block format as (P*xt).
  void permutedFFT(const N_LAS_BlockVector & xt, N_LAS_BlockVector * xf); // Reference for input, pointer for output
  // xt = P^{-1}D^{-1}*xf
  void permutedIFT(const N_LAS_BlockVector & xf, N_LAS_BlockVector * xt);

  // Registration method for the device packaage
  void registerAppLoader( Teuchos::RCP<N_LOA_Loader> appLoaderPtr );
  
  // Register the interface to the devices
  void registerDeviceInterface( Teuchos::RCP< N_DEV_DeviceInterface > devIntPtr );

  // Registration method for the Application Vector
  void registerAppVec( Teuchos::RCP<N_LAS_Vector> appVecPtr );

  // Registration method for the next Application State Vector
  void registerAppNextStaVec( Teuchos::RCP<N_LAS_Vector> appNextStaVecPtr );

  // Registration method for the curr Application State Vector
  void registerAppCurrStaVec( Teuchos::RCP<N_LAS_Vector> appCurrStaVecPtr );

  // Registration method for the last Application State Vector
  void registerAppLastStaVec( Teuchos::RCP<N_LAS_Vector> appLastStaVecPtr );

  // Registration method for the next Application Store Vector
  void registerAppNextStoVec( Teuchos::RCP<N_LAS_Vector> appNextStoVecPtr );

  // Registration method for the curr Application Store Vector
  void registerAppCurrStoVec( Teuchos::RCP<N_LAS_Vector> appCurrStoVecPtr );

  // Registration method for the last Application Store Vector
  void registerAppLastStoVec( Teuchos::RCP<N_LAS_Vector> appLastStoVecPtr );
  
  // Registration method for the Application Store Vector Lead Current Q Component 
  void registerAppStoLeadCurrQCompVec( RefCountPtr<N_LAS_Vector> appStoLeadCurrQCompVecPtr );

  // Registration method for the Application dQdx Matrix
  void registerAppdQdx( Teuchos::RCP<N_LAS_Matrix> appdQdxPtr );

  // Registration method for the Application dFdx Matrix
  void registerAppdFdx( Teuchos::RCP<N_LAS_Matrix> appdFdxPtr );

  // Registration method for the MPDE size dQdx Block Matrix
  void registerMPDEdQdx( Teuchos::RCP<N_LAS_BlockMatrix> bmdQdxPtr );

  // Registration method for the MPDE size dFdx Block Matrix
  void registerMPDEdFdx( Teuchos::RCP<N_LAS_BlockMatrix> bmdFdxPtr );

  // Registration method for the HB size dQdx Block Matrix
  void registerHBdQdx( Teuchos::RCP<N_LAS_BlockMatrix> bmdQdxPtr );

  // Registration method for the HB size dFdx Block Matrix
  void registerHBdFdx( Teuchos::RCP<N_LAS_BlockMatrix> bmdFdxPtr );

  void registerXt(Teuchos::RCP<N_LAS_BlockVector> XtPtr);

  void registerVt(Teuchos::RCP<N_LAS_BlockVector> VtPtr);
  
  void registerHBBuilder(Teuchos::RCP<N_LAS_HBBuilder> hbBuilder)
  { hbBuilderRCPtr_ = hbBuilder; }

  private :

  //MPDE State
  N_MPDE_State & state_;

  // discretization
  Teuchos::RCP<const N_MPDE_Discretization> fastTimeDiscRCPtr_;

  //Fast Time Scale Points
  vector<double> times_;
  int periodicTimesOffset_;
  vector<double> periodicTimes_;
  double period_;
 
  // Matrix free flag, operator is being applied not loaded
  bool matrixFreeFlag_;
 
  //a vector of bools to indicate if a solution variable does not
  //appear to be periodic.  Non-periodic signals will be handled
  //differently by the MPDE_Loader class
  vector<bool> nonPeriodic_;

  // Base Application loader
  Teuchos::RCP<N_LOA_Loader> appLoaderPtr_; // Actually a CktLoader
  Teuchos::RCP< N_DEV_DeviceInterface > devInterfacePtr_;

  // Application Linear Objects
  Teuchos::RCP<N_LAS_Vector> appVecPtr_;
  Teuchos::RCP<N_LAS_Vector> appNextStaVecPtr_;
  Teuchos::RCP<N_LAS_Vector> appCurrStaVecPtr_;
  Teuchos::RCP<N_LAS_Vector> appLastStaVecPtr_;
  Teuchos::RCP<N_LAS_Matrix> appdQdxPtr_;
  Teuchos::RCP<N_LAS_Matrix> appdFdxPtr_;

  Teuchos::RCP<N_LAS_Vector> appNextStoVecPtr_;
  Teuchos::RCP<N_LAS_Vector> appCurrStoVecPtr_;
  Teuchos::RCP<N_LAS_Vector> appLastStoVecPtr_;
  Teuchos::RCP<N_LAS_Vector> appStoLeadCurrQCompVecPtr_;

  // Tmp storage block matrices 
  Teuchos::RCP<N_LAS_BlockMatrix> bmdQdxPtr_;
  Teuchos::RCP<N_LAS_BlockMatrix> bmdFdxPtr_;

  // Additional Tmp storage block matrices  (HB case)
  Teuchos::RCP<N_LAS_BlockMatrix> bmdQdx2Ptr_;
  Teuchos::RCP<N_LAS_BlockMatrix> bmdFdx2Ptr_;

  // HB Builder:  (needed to convert AztecOO created N_LAS_Vectors into N_LAS_BlockVectors
  Teuchos::RCP<N_LAS_HBBuilder> hbBuilderRCPtr_;

//  Ref CountPtr<N_LAS_BlockVector> bXtPtr_;
  Teuchos::RCP<N_LAS_BlockVector> bXtPtr_;
  Teuchos::RCP<N_LAS_BlockVector> bVtPtr_;   

};


//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerAppLoader
// Purpose       : Registration method for the device packaage
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerAppLoader( Teuchos::RCP<N_LOA_Loader> appLoaderPtr )
{ 
  appLoaderPtr_ = appLoaderPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerDeviceInterface
// Purpose       : Registration method for the device packaage
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerDeviceInterface( Teuchos::RCP< N_DEV_DeviceInterface > devInterfacePtr )
{ 
  devInterfacePtr_ = devInterfacePtr;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerAppVec
// Purpose       : Registration method for the Application Vector
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerAppVec( Teuchos::RCP<N_LAS_Vector> appVecPtr )
{ 
  appVecPtr_ = appVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerAppNextStaVec
// Purpose       : Registration method for the Application State Vector
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerAppNextStaVec( Teuchos::RCP<N_LAS_Vector> appNextStaVecPtr )
{ 
  appNextStaVecPtr_ = appNextStaVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerAppCurrStaVec
// Purpose       : Registration method for the Application State Vector
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 01/25/07
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerAppCurrStaVec( Teuchos::RCP<N_LAS_Vector> appCurrStaVecPtr )
{ 
  appCurrStaVecPtr_ = appCurrStaVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerAppLastStaVec
// Purpose       : Registration method for the Application State Vector
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 01/25/07
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerAppLastStaVec( Teuchos::RCP<N_LAS_Vector> appLastStaVecPtr )
{ 
  appLastStaVecPtr_ = appLastStaVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerAppNextStoVec
// Purpose       : Registration method for the Application Store Vector
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerAppNextStoVec( Teuchos::RCP<N_LAS_Vector> appNextStoVecPtr )
{ 
  appNextStoVecPtr_ = appNextStoVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerAppCurrStoVec
// Purpose       : Registration method for the Application Store Vector
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 01/25/07
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerAppCurrStoVec( Teuchos::RCP<N_LAS_Vector> appCurrStoVecPtr )
{ 
  appCurrStoVecPtr_ = appCurrStoVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerAppLastStoVec
// Purpose       : Registration method for the Application Store Vector
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 01/25/07
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerAppLastStoVec( Teuchos::RCP<N_LAS_Vector> appLastStoVecPtr )
{ 
  appLastStoVecPtr_ = appLastStoVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerAppStoLeadCurrQCompVec
// Purpose       : Registration method for the Application Store Vector
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek
// Creation Date : 
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerAppStoLeadCurrQCompVec( RefCountPtr<N_LAS_Vector> appStoLeadCurrQCompVecPtr )
{ 
  appStoLeadCurrQCompVecPtr_ = appStoLeadCurrQCompVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerXt
// Purpose       : Registration method for the Application dQdx Matrix
// Special Notes :
// Scope         : public
// Creator       : Ting  Mei
// Creation Date : 09/04/08
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerXt(Teuchos::RCP<N_LAS_BlockVector> XtPtr)
{ 
  bXtPtr_ = XtPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerVt
// Purpose       : Registration method for the Application dQdx Matrix
// Special Notes :
// Scope         : public
// Creator       : Ting  Mei
// Creation Date : 09/04/08
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerVt(Teuchos::RCP<N_LAS_BlockVector> VtPtr)
{ 
  bVtPtr_ = VtPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerAppdQdx
// Purpose       : Registration method for the Application dQdx Matrix
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerAppdQdx( Teuchos::RCP<N_LAS_Matrix> appdQdxPtr )
{ 
  appdQdxPtr_ = appdQdxPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerAppdFdx
// Purpose       :  Registration method for the Application dFdx Matrix
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_LOA_HBLoader::registerAppdFdx( Teuchos::RCP<N_LAS_Matrix> appdFdxPtr )
{ 
  appdFdxPtr_ = appdFdxPtr; 
}

// Registration method for the MPDE size dQdx Block Matrix
inline void N_LOA_HBLoader::registerMPDEdQdx( Teuchos::RCP<N_LAS_BlockMatrix> bmdQdxPtr )
{ 
  bmdQdxPtr_ = bmdQdxPtr; 
}

// Registration method for the MPDE size dFdx Block Matrix
inline void N_LOA_HBLoader::registerMPDEdFdx( Teuchos::RCP<N_LAS_BlockMatrix> bmdFdxPtr )
{ 
  bmdFdxPtr_ = bmdFdxPtr; 
}

// Registration method for the HB size dQdx Block Matrix
inline void N_LOA_HBLoader::registerHBdQdx( Teuchos::RCP<N_LAS_BlockMatrix> bmdQdxPtr )
{ 
  bmdQdx2Ptr_ = bmdQdxPtr; 
}

// Registration method for the HB size dFdx Block Matrix
inline void N_LOA_HBLoader::registerHBdFdx( Teuchos::RCP<N_LAS_BlockMatrix> bmdFdxPtr )
{ 
  bmdFdx2Ptr_ = bmdFdxPtr; 
}

#endif // Xyce_LOA_HBLoader_H


