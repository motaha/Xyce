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
// Filename       : $RCSfile: N_LAS_System.h,v $
//
// Purpose        : Container class for linear system Jacobian, RHS, Soln Vec,
//                  Error Vec, and Creators N_LAS_QueryUtil and N_PDS_ParMap
//                  are registered to minimize input for matrix and vector
//                  creation
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/10/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.57 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  Xyce_LAS_SYSTEM_H
#define  Xyce_LAS_SYSTEM_H

// ---------- Standard Includes ----------

#include <list>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>
#include <N_ANP_fwd.h>
// to eliminate RCP warnings, putting N_PDS_Manager header here.
#include <N_PDS_Manager.h>

// ---------- Forward Declarations ----------

class N_LAS_QueryUtil;
class N_LAS_LAFactory;
class N_LAS_Builder;
class N_LAS_Vector;
class N_LAS_MultiVector;
class N_LAS_Matrix;
class N_LAS_IterativeSolver;

class Epetra_MapColoring;

//-----------------------------------------------------------------------------
// Class         : N_LAS_System
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
class N_LAS_System
{

public:

  // Default Constructor
  N_LAS_System()
  : pdsMgr_(0),
    lasQueryUtil_(0),
    lasBuilder_(0),
    anpIntPtr_(0),
    jacobianMatrixPtr_(0),
    rhsVectorPtr_(0),
    jdxpVectorPtr_(0),
    fVectorPtr_(0),

    daeQVectorPtr_(0),
    daeFVectorPtr_(0),
    daeBVectorPtr_(0),
    dFdxdVpVectorPtr_(0),
    dQdxdVpVectorPtr_(0),
    dQdxMatrixPtr_(0),
    dFdxMatrixPtr_(0),
    currSolVectorPtrPtr_(0),
    nextSolVectorPtrPtr_(0),
    lastSolVectorPtrPtr_(0),
    flagSolVectorPtrPtr_(0),
    nextSolDerivVectorPtrPtr_(0),
    nextStaDerivVectorPtrPtr_(0),
    currStaVectorPtrPtr_(0),
    nextStaVectorPtrPtr_(0),
    lastStaVectorPtrPtr_(0),
    currStoVectorPtrPtr_(0),
    nextStoVectorPtrPtr_(0),
    lastStoVectorPtrPtr_(0),
    tmpStaVectorPtrPtr_(0),
    tmpSolVectorPtrPtr_(0),
    tmpStaDerivVectorPtrPtr_(0),
    tmpSolDerivVectorPtrPtr_(0),
    tmpStaDivDiffVectorPtrPtr_(0),
    tmpSolDivDiffVectorPtrPtr_(0) ,
    deviceMaskVectorPtrPtr_(0) ,
    solnColoringPtr_(0),
    initialConditionColoringPtr_(0),
    nonTrivialDeviceMaskFlag_(false)
  {}
  
  //Destructor
  ~N_LAS_System();

  // Registration methods for necessary utilities

  // Registers the Parallel Distribution Services (PDS) manager (i.e., sets the
  // pointer.
  bool registerPDSManager(N_PDS_Manager * PDS_Manager)
  { return ( pdsMgr_ = PDS_Manager ); }

  bool registerANPInterface(N_ANP_AnalysisInterface * ANP_Interface)
  { return ( anpIntPtr_ = ANP_Interface ); }

  // Registers the Query utility
  bool registerQueryUtil(N_LAS_QueryUtil * LAS_QUtil)
  { return ( lasQueryUtil_ = LAS_QUtil ); }

  // Registers the LAS Builder object
  bool registerBuilder(N_LAS_Builder * builder)
  { return ( lasBuilder_ = builder ); }

  // Registration methods for soln and state vector objects since the TIA or
  // NLS packages will probably own these

  // Registers the pointer to the nonlinear system's Jacobian matrix
  bool registerJacobianMatrix(N_LAS_Matrix * jacMatrixPtr)
  { return ( jacobianMatrixPtr_ = jacMatrixPtr ); }

#ifdef Xyce_DEBUG_VOLTLIM
  // Registers the pointer to the Jacobian test matrix (old DAE)
  bool registerJacTestMatrix(N_LAS_Matrix * jacTestMatrixPtr)
  { return ( jacTestMatrixPtr_ = jacTestMatrixPtr ); }

  // Registers the pointer to the Jacobian test matrix (new DAE)
  bool registerdFdxTestMatrix(N_LAS_Matrix * dFdxTestMatrixPtr)
  { return ( dFdxTestMatrixPtr_ = dFdxTestMatrixPtr ); }
  bool registerdQdxTestMatrix(N_LAS_Matrix * dQdxTestMatrixPtr)
  { return ( dQdxTestMatrixPtr_ = dQdxTestMatrixPtr ); }
#endif

  // Registers the pointer to the nonlinear system's residual (RHS) vector
  bool registerRHSVector (N_LAS_Vector * rhsVecPtr)
  { return ( rhsVectorPtr_ = rhsVecPtr ); }

  // Registers the pointer to the nonlinear system's limiting (jdxp) vector
  bool registerJDXPVector (N_LAS_Vector * jdxpVecPtr)
  { return ( jdxpVectorPtr_ = jdxpVecPtr ); }

  // Registers the pointer to the nonlinear system's F vector
  bool registerFVector   (N_LAS_Vector * fVecPtr)
  { return ( fVectorPtr_ = fVecPtr ); }

#ifdef Xyce_DEBUG_VOLTLIM
  // Registers the pointer to the voltlim dx vector
  bool registerDxVoltlimVector   (N_LAS_Vector * dxVLPtr)
  { return ( dxVoltlimVectorPtr_ = dxVLPtr ); }

  // Registers the pointer to the nonlinear system's limiting (jdx2) vector (old DAE)
  bool registerJDX2Vector (N_LAS_Vector * jdx2VecPtr)
  { return ( jdx2VectorPtr_ = jdx2VecPtr ); }

  // Registers the pointer to the nonlinear system's limiting vectors (new DAE)
  bool registerFDX2Vector (N_LAS_Vector * Fdx2VecPtr)
  { return ( Fdx2VectorPtr_ = Fdx2VecPtr ); }

  bool registerQDX2Vector (N_LAS_Vector * Qdx2VecPtr)
  { return ( Qdx2VectorPtr_ = Qdx2VecPtr ); }
#endif
  // DAE vector registrations:
  bool registerDAEQVector (N_LAS_Vector * qVecPtr)
  {  return (daeQVectorPtr_ = qVecPtr); }

  bool registerDAEFVector (N_LAS_Vector * fVecPtr)
  {  return (daeFVectorPtr_ = fVecPtr); }

  bool registerDAEBVector (N_LAS_Vector * bVecPtr)
  {  return (daeBVectorPtr_ = bVecPtr); }

  // dFdxdvp vector, for new-DAE voltage limiting.
  bool registerdFdxdVpVector (N_LAS_Vector * dFdxdVpVectorPtr)
  { return ( dFdxdVpVectorPtr_ = dFdxdVpVectorPtr); }

  // dQdxdvp vector, for new-DAE voltage limiting.
  bool registerdQdxdVpVector (N_LAS_Vector * dQdxdVpVectorPtr)
  { return ( dQdxdVpVectorPtr_ = dQdxdVpVectorPtr) ; }

  // DAE matrix registrations:
  bool registerDAEdQdxMatrix (N_LAS_Matrix * dqdxMatPtr)
  {  return (dQdxMatrixPtr_ = dqdxMatPtr); }

  bool registerDAEdFdxMatrix (N_LAS_Matrix * dfdxMatPtr)
  {  return (dFdxMatrixPtr_ = dfdxMatPtr); }

  // Registers the pointer to the current solution vector
  bool registerCurrSolVector(N_LAS_Vector ** solVecPtrPtr)
  { return ( currSolVectorPtrPtr_ = solVecPtrPtr ); }

  // Registers the pointer to the next solution vector
  bool registerNextSolVector(N_LAS_Vector ** solVecPtrPtr)
  { return ( nextSolVectorPtrPtr_ = solVecPtrPtr ); }

  // Registers the pointer to the flag solution vector
  bool registerFlagSolVector(N_LAS_Vector ** flagSolVecPtrPtr)
  { return ( flagSolVectorPtrPtr_ = flagSolVecPtrPtr ); }

  // Registers the pointer to the previous solution vector
  bool registerLastSolVector(N_LAS_Vector ** solVecPtrPtr)
  { return ( lastSolVectorPtrPtr_ = solVecPtrPtr ); }

  // Registers the pointer to the next solution vector's derivatives
  bool registerNextSolDerivVector(N_LAS_Vector ** solVecPtrPtr)
  { return ( nextSolDerivVectorPtrPtr_ = solVecPtrPtr ); }

  // Registers the pointer to the next state (auxiliary) vector's derivatives
  bool registerNextStaDerivVector(N_LAS_Vector ** staVecPtrPtr)
  { return ( nextStaDerivVectorPtrPtr_ = staVecPtrPtr ); }

  // Registers the pointer to the current state (auxiliary) vector
  bool registerCurrStaVector(N_LAS_Vector ** stateVecPtrPtr)
  { return ( currStaVectorPtrPtr_ = stateVecPtrPtr ); }

  // Registers the pointer to the next state (auxiliary) vector
  bool registerNextStaVector(N_LAS_Vector ** stateVecPtrPtr)
  { return ( nextStaVectorPtrPtr_ = stateVecPtrPtr ); }

  // Registers the pointer to the previous state (auxiliary) vector
  bool registerLastStaVector(N_LAS_Vector ** stateVecPtrPtr)
  { return ( lastStaVectorPtrPtr_ = stateVecPtrPtr ); }

  // Registers the pointer to the current store (auxiliary) vector
  bool registerCurrStoVector(N_LAS_Vector ** storeVecPtrPtr)
  { return ( currStoVectorPtrPtr_ = storeVecPtrPtr ); }

  // Registers the pointer to the next store (auxiliary) vector
  bool registerNextStoVector(N_LAS_Vector ** storeVecPtrPtr)
  { return ( nextStoVectorPtrPtr_ = storeVecPtrPtr ); }

  // Registers the pointer to the previous store (auxiliary) vector
  bool registerLastStoVector(N_LAS_Vector ** storeVecPtrPtr)
  { return ( lastStoVectorPtrPtr_ = storeVecPtrPtr ); }

  // Registers the pointer to the temporary state (auxiliary) vector
  bool registerTmpStaVector(N_LAS_Vector ** tmpStaVecPtrPtr)
  { return ( tmpStaVectorPtrPtr_ = tmpStaVecPtrPtr ); }

  // Registers the pointer to the temporary solution vector
  bool registerTmpSolVector(N_LAS_Vector ** tmpSolVecPtrPtr)
  { return ( tmpSolVectorPtrPtr_ = tmpSolVecPtrPtr ); }

  // Registers the pointer to the temporary state (auxiliary) vector's
  // derivatives
  bool registerTmpStaDerivVector(N_LAS_Vector ** tmpStaDerivVecPtrPtr)
  { return ( tmpStaDerivVectorPtrPtr_ = tmpStaDerivVecPtrPtr ); }

  // Registers the pointer to the temporary solution vector's derivatives
  bool registerTmpSolDerivVector(N_LAS_Vector ** tmpSolDerivVecPtrPtr)
  { return ( tmpSolDerivVectorPtrPtr_ = tmpSolDerivVecPtrPtr ); }

  bool registerTmpStaDivDiffVector(N_LAS_Vector ** tmpStaDivDiffVecPtrPtr)
  { return ( tmpSolDivDiffVectorPtrPtr_ = tmpStaDivDiffVecPtrPtr ); }
  bool registerTmpSolDivDiffVector(N_LAS_Vector ** tmpSolDivDiffVecPtrPtr)
  { return ( tmpSolDivDiffVectorPtrPtr_ = tmpSolDivDiffVecPtrPtr ); }

  // Register pointer to device mask
  bool registerDeviceMaskVector(N_LAS_Vector ** dmPP)
  { return (deviceMaskVectorPtrPtr_ = dmPP); }
  
  bool registerSolnColoring( Epetra_MapColoring * coloring )
  { return ( solnColoringPtr_ = coloring ); }

  bool registerICColoring( Epetra_MapColoring * coloring )
  { return ( initialConditionColoringPtr_ = coloring ); }

  // set boolean telling whether the device mask has any zeros  
  void setNonTrivialDeviceMaskFlag(bool nTDMF)
  { nonTrivialDeviceMaskFlag_ = nTDMF; }

  // Accessor methods for linear system objects

  // Get method for the Jacobian matrix
  N_LAS_Matrix *  getJacobianMatrix() { return jacobianMatrixPtr_; }

#ifdef Xyce_DEBUG_VOLTLIM
  // Get method for the Jacobian test matrix (old DAE)
  N_LAS_Matrix *  getJacTestMatrix() { return jacTestMatrixPtr_; }
  // Get method for the Jacobian test matrix (new DAE)
  N_LAS_Matrix *  getdFdxTestMatrix() { return dFdxTestMatrixPtr_; }
  N_LAS_Matrix *  getdQdxTestMatrix() { return dQdxTestMatrixPtr_; }
#endif

  // Get method for the residual (RHS) vector
  N_LAS_Vector *  getRHSVector() { return rhsVectorPtr_; }

  // Get method for the limiting (jdxp) vector
  N_LAS_Vector *  getJDXPVector() { return jdxpVectorPtr_; }

  // Get method for the F vector
  N_LAS_Vector *  getFVector() { return fVectorPtr_; }

#ifdef Xyce_DEBUG_VOLTLIM
  // Get method for the voltage limiting dx vector
  N_LAS_Vector *  getDxVoltlimVector() { return dxVoltlimVectorPtr_; }

  // Get method for the limiting (jdx2) vector (old DAE)
  N_LAS_Vector *  getJDX2Vector() { return jdx2VectorPtr_; }

  // Get methods for the limiting (jdx2) vector (new DAE)
  N_LAS_Vector *  getFDX2Vector() { return Fdx2VectorPtr_; }
  N_LAS_Vector *  getQDX2Vector() { return Qdx2VectorPtr_; }
#endif

  // Get methods for DAE vectors:
  N_LAS_Vector * getDAEQVector () { return daeQVectorPtr_; }

  N_LAS_Vector * getDAEFVector () { return daeFVectorPtr_; }

  N_LAS_Vector * getDAEBVector () { return daeBVectorPtr_; }

  // dFdxdvp vector, for new-DAE voltage limiting.
  N_LAS_Vector * getdFdxdVpVector () { return  dFdxdVpVectorPtr_; }

  // dQdxdvp vector, for new-DAE voltage limiting.
  N_LAS_Vector * getdQdxdVpVector () { return  dQdxdVpVectorPtr_; }

  // Get methods for DAE matrices:
  N_LAS_Matrix * getDAEdQdxMatrix () { return dQdxMatrixPtr_; }

  N_LAS_Matrix * getDAEdFdxMatrix () { return dFdxMatrixPtr_; }

  // Get method for the current solution vector
  N_LAS_Vector *  getCurrSolVector() { return *currSolVectorPtrPtr_; }

  // Get method for the current solution vector pointer
  N_LAS_Vector ** getCurrSolVectorPtr() { return currSolVectorPtrPtr_; }

  // Get method for next current solution vector
  N_LAS_Vector *  getNextSolVector() { return *nextSolVectorPtrPtr_; }

  // Get method for next current solution vector pointer
  N_LAS_Vector ** getNextSolVectorPtr() { return nextSolVectorPtrPtr_; }

  // Get method for previous current solution vector 
  N_LAS_Vector *  getLastSolVector() { return *lastSolVectorPtrPtr_; }

  // Get method for previous current solution vector pointer
  N_LAS_Vector ** getLastSolVectorPtr() { return lastSolVectorPtrPtr_; }

  // Get method for flag solution vector
  N_LAS_Vector * getFlagSolVector() { return *flagSolVectorPtrPtr_; }

  // Get method for flag solution vector
  N_LAS_Vector ** getFlagSolVectorPtr() { return flagSolVectorPtrPtr_; }

  N_LAS_Vector *  getNextSolDerivVector()
  {
    return *nextSolDerivVectorPtrPtr_;
  }
  N_LAS_Vector ** getNextSolDerivVectorPtr()
  {
    return nextSolDerivVectorPtrPtr_;
  }

  N_LAS_Vector *  getCurrStaVector() { return *currStaVectorPtrPtr_; }
  N_LAS_Vector ** getCurrStaVectorPtr() { return currStaVectorPtrPtr_; }
  N_LAS_Vector *  getNextStaVector() { return *nextStaVectorPtrPtr_; }
  N_LAS_Vector ** getNextStaVectorPtr() { return nextStaVectorPtrPtr_; }
  N_LAS_Vector *  getLastStaVector () { return *lastStaVectorPtrPtr_; }
  N_LAS_Vector ** getLastStaVectorPtr() { return lastStaVectorPtrPtr_; }
  N_LAS_Vector *  getCurrStoVector() { return *currStoVectorPtrPtr_; }
  N_LAS_Vector ** getCurrStoVectorPtr() { return currStoVectorPtrPtr_; }
  N_LAS_Vector *  getNextStoVector() { return *nextStoVectorPtrPtr_; }
  N_LAS_Vector ** getNextStoVectorPtr() { return nextStoVectorPtrPtr_; }
  N_LAS_Vector *  getLastStoVector () { return *lastStoVectorPtrPtr_; }
  N_LAS_Vector ** getLastStoVectorPtr() { return lastStoVectorPtrPtr_; }

  N_LAS_Vector *  getNextStaDerivVector()
  {
    return *nextStaDerivVectorPtrPtr_;
  }
  N_LAS_Vector ** getNextStaDerivVectorPtr()
  {
    return nextStaDerivVectorPtrPtr_;
  }

  N_LAS_Vector *  getTmpSolVector   () { return *tmpSolVectorPtrPtr_; }
  N_LAS_Vector ** getTmpSolVectorPtr() { return tmpSolVectorPtrPtr_; }
  N_LAS_Vector *  getTmpStaVector   () { return *tmpStaVectorPtrPtr_; }
  N_LAS_Vector ** getTmpStaVectorPtr() { return tmpStaVectorPtrPtr_; }

  N_LAS_Vector ** getTmpSolDerivVectorPtr()
  {
    return tmpSolDerivVectorPtrPtr_;
  }
  N_LAS_Vector ** getTmpStaDerivVectorPtr()
  {
    return tmpStaDerivVectorPtrPtr_;
  }

  N_LAS_Vector ** getTmpSolDivDiffVectorPtr()
  {
    return tmpSolDivDiffVectorPtrPtr_;
  }
  N_LAS_Vector ** getTmpStaDivDiffVectorPtr()
  {
    return tmpStaDivDiffVectorPtrPtr_;
  }

  // Get method for the device mask vector
  N_LAS_Vector *  getDeviceMaskVector() { return *deviceMaskVectorPtrPtr_; }
  N_LAS_Vector **  getDeviceMaskVectorPtr() { return deviceMaskVectorPtrPtr_; }
  // get boolean telling whether the device mask has any zeros  
  bool getNonTrivialDeviceMaskFlag() {return nonTrivialDeviceMaskFlag_; }

  // Registers the Query utility
  N_LAS_QueryUtil * getQueryUtil()
  { 
    return lasQueryUtil_; 
  }
    
  // Returns the Parallel Distribution Services (PDS) manager 
  N_PDS_Manager *  getPDSManager()
  { 
    return pdsMgr_; 
  }

  // Problem size accessors

  // Get method for the GLOBAL solution vector size (i.e., number of unknowns)
  int getGlobalSolutionSize();

  // Get method for the GLOBAL state vector size (i.e., number of unknowns)
  int getGlobalStateSize();

  // Get method for the solution vector size (i.e., number of unknowns)
  int getSolutionSize();

  // Get method for the state (auxiliary) vector size
  int getStateSize();

  // Get method for the residual (RHS) vector size
  int getRHSSize();

  // Create residual (RHS) and Jacobian
  bool initializeSystem();

  // Builder access
  N_LAS_Builder & builder() { return *lasBuilder_; }

  bool updateExternValsSolnVector(N_LAS_MultiVector * solnVector);
  bool updateExternValsStateVector(N_LAS_MultiVector * stateVector);
  bool updateExternValsStoreVector(N_LAS_MultiVector * storeVector);

  bool localInMatrix(const int & row) const;
  bool localInSolnVector(const int & row) const;
  bool localInStateVector(const int & row) const;

  // Debug call outputs vectors to console
  void debug() const;

private:

  N_PDS_Manager   *  pdsMgr_;

  N_LAS_QueryUtil *  lasQueryUtil_;
  N_LAS_Builder   *  lasBuilder_;

  N_ANP_AnalysisInterface * anpIntPtr_;

  N_LAS_Matrix    *  jacobianMatrixPtr_;
  N_LAS_Vector    *  rhsVectorPtr_;
  N_LAS_Vector    *  jdxpVectorPtr_;
  N_LAS_Vector    *  fVectorPtr_;
#ifdef Xyce_DEBUG_VOLTLIM
  N_LAS_Vector    *  dxVoltlimVectorPtr_;

  // old dae:
  N_LAS_Vector    *  jdx2VectorPtr_;
  N_LAS_Matrix    *  jacTestMatrixPtr_;

  // new dae:
  N_LAS_Vector    *  Fdx2VectorPtr_;
  N_LAS_Matrix    *  dFdxTestMatrixPtr_;
  N_LAS_Vector    *  Qdx2VectorPtr_;
  N_LAS_Matrix    *  dQdxTestMatrixPtr_;
#endif
  // DAE vectors:
  N_LAS_Vector    *  daeQVectorPtr_;
  N_LAS_Vector    *  daeFVectorPtr_;
  N_LAS_Vector    *  daeBVectorPtr_;

  N_LAS_Vector    *  dFdxdVpVectorPtr_;
  N_LAS_Vector    *  dQdxdVpVectorPtr_;

  // DAE matrices:
  N_LAS_Matrix    * dQdxMatrixPtr_;
  N_LAS_Matrix    * dFdxMatrixPtr_;

  N_LAS_Vector ** currSolVectorPtrPtr_;
  N_LAS_Vector ** nextSolVectorPtrPtr_;
  N_LAS_Vector ** lastSolVectorPtrPtr_;
  N_LAS_Vector ** flagSolVectorPtrPtr_;

  N_LAS_Vector ** nextSolDerivVectorPtrPtr_;
  N_LAS_Vector ** nextStaDerivVectorPtrPtr_;

  N_LAS_Vector ** currStaVectorPtrPtr_;
  N_LAS_Vector ** nextStaVectorPtrPtr_;
  N_LAS_Vector ** lastStaVectorPtrPtr_;
  N_LAS_Vector ** currStoVectorPtrPtr_;
  N_LAS_Vector ** nextStoVectorPtrPtr_;
  N_LAS_Vector ** lastStoVectorPtrPtr_;

  N_LAS_Vector **  tmpStaVectorPtrPtr_;
  N_LAS_Vector **  tmpSolVectorPtrPtr_;

  N_LAS_Vector **  tmpStaDerivVectorPtrPtr_;
  N_LAS_Vector **  tmpSolDerivVectorPtrPtr_;

  N_LAS_Vector **  tmpStaDivDiffVectorPtrPtr_;
  N_LAS_Vector **  tmpSolDivDiffVectorPtrPtr_;

  N_LAS_Vector **  deviceMaskVectorPtrPtr_;

  Epetra_MapColoring * solnColoringPtr_;

  Epetra_MapColoring * initialConditionColoringPtr_;

  bool nonTrivialDeviceMaskFlag_;
};

#endif

