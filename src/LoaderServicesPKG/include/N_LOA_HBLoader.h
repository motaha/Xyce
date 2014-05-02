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
// Revision Number: $Revision: 1.10 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
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
class N_LAS_Builder;

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
                    const Teuchos::RCP<const N_MPDE_Discretization> discPtr
                  )
  : state_(state),
    fastTimeDiscPtr_(discPtr),
    periodicTimesOffset_(0),
    period_(0),
    matrixFreeFlag_(false)
  {}

  // Destructor
  ~N_LOA_HBLoader() {}

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

  // Get the voltage limiter flag:
  bool getLimiterFlag () { return N_LOA_HBLoader::appLoaderPtr_->getLimiterFlag (); }

  // Get the stored time-domain Jacobians from the HB loader.
  std::vector<Teuchos::RCP<N_LAS_Matrix> >& getStoredQdx() { return vecAppdQdxPtr_; }
  std::vector<Teuchos::RCP<N_LAS_Matrix> >& getStoredFdx() { return vecAppdFdxPtr_; }

  // New functions for LOA_HBLoader
  // Assumption:  xt is in the block format:
  // xt.block(i) = { x_n(t_i) }_{n=0..N} N = number of solution components, t_i = time point i
  // P takes xt's block format and converts it to: 
  // (P*xt).block(n) = { x_n(t_i) }_{i=0..T} T = number of time points in block vector.

   
  void setFastTimes( const std::vector<double> & times );
  void constructPeriodicTimes();

  Teuchos::RCP<N_LAS_BlockVector> & getStoreVecFreqPtr()  { return bStoreVecFreqPtr_;} 
  // Is the problem being applied as a matrix free operator.
  bool matrixFreeFlag() const { return matrixFreeFlag_; }
  void setMatrixFreeFlag(bool matrixFreeFlag) { matrixFreeFlag_ = matrixFreeFlag; }
 
  // xf = D*P*xt, xf has the same block format as (P*xt).
  void permutedFFT(const N_LAS_BlockVector & xt, N_LAS_BlockVector * xf); // Reference for input, pointer for output
  // xt = P^{-1}D^{-1}*xf
  void permutedIFT(const N_LAS_BlockVector & xf, N_LAS_BlockVector * xt);

  // Registration method for the device packaage
  void registerAppLoader( Teuchos::RCP<N_LOA_Loader> appLoaderPtr )
  { appLoaderPtr_ = appLoaderPtr; }
  
  // Register the interface to the devices
  void registerDeviceInterface( Teuchos::RCP< N_DEV_DeviceInterface > devIntPtr )
  { devInterfacePtr_ = devIntPtr; }

  void registerHBBuilder(Teuchos::RCP<N_LAS_HBBuilder> hbBuilderPtr);

  void registerAppBuilder(Teuchos::RCP<N_LAS_Builder> appBuilderPtr);

    virtual bool setParam (std::string & name, double val) { return false; }
    virtual double getParamAndReduce (const std::string & name) { return 0.0; }
    virtual bool getParamAndReduce (const std::string & name, double & val) { return false; }
    
  private :

  //MPDE State
  N_MPDE_State & state_;

  // discretization
  Teuchos::RCP<const N_MPDE_Discretization> fastTimeDiscPtr_;

  //Fast Time Scale Points
  std::vector<double> times_;
  int periodicTimesOffset_;
  std::vector<double> periodicTimes_;
  double period_;
 
  // Matrix free flag, operator is being applied not loaded
  bool matrixFreeFlag_;
 
  // Base Application loader
  Teuchos::RCP<N_LOA_Loader> appLoaderPtr_; // Actually a CktLoader
  Teuchos::RCP< N_DEV_DeviceInterface > devInterfacePtr_;

  // Application Linear Objects
  Teuchos::RCP<N_LAS_Vector> appVecPtr_;
  Teuchos::RCP<N_LAS_Vector> appNextStaVecPtr_;
  Teuchos::RCP<N_LAS_Vector> appCurrStaVecPtr_;
  Teuchos::RCP<N_LAS_Vector> appLastStaVecPtr_;
  Teuchos::RCP<N_LAS_Matrix> appdQdxPtr_;
  std::vector<Teuchos::RCP<N_LAS_Matrix> > vecAppdQdxPtr_;
  Teuchos::RCP<N_LAS_Matrix> appdFdxPtr_;
  std::vector<Teuchos::RCP<N_LAS_Matrix> > vecAppdFdxPtr_;

  Teuchos::RCP<N_LAS_Vector> appNextStoVecPtr_;
  Teuchos::RCP<N_LAS_Vector> appCurrStoVecPtr_;
  Teuchos::RCP<N_LAS_Vector> appLastStoVecPtr_;
  Teuchos::RCP<N_LAS_Vector> appStoLeadCurrQCompVecPtr_;

  // Tmp storage block matrices 
  Teuchos::RCP<N_LAS_BlockMatrix> bmdQdxPtr_;
  Teuchos::RCP<N_LAS_BlockMatrix> bmdFdxPtr_;

  // HB Builder:  (needed to convert AztecOO created N_LAS_Vectors into N_LAS_BlockVectors
  Teuchos::RCP<N_LAS_HBBuilder> hbBuilderPtr_;

  // App Builder:  (needed to load time domain vectors and matrices)
  Teuchos::RCP<N_LAS_Builder> appBuilderPtr_;

  Teuchos::RCP<N_LAS_BlockVector> bXtPtr_;
  Teuchos::RCP<N_LAS_BlockVector> bVtPtr_;

  Teuchos::RCP<N_LAS_BlockVector> bStoreVecFreqPtr_;
  Teuchos::RCP<N_LAS_BlockVector> bStoreLeadCurrQCompVecFreqPtr_; 

};

#endif // Xyce_LOA_HBLoader_H


