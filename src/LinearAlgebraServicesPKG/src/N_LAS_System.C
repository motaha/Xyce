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
// Filename       : $RCSfile: N_LAS_System.C,v $
//
// Purpose        : Container class for linear system
//                  Jacobian, RHS, Soln Vec, Error Vec, and Creators
//                  N_LAS_QueryUtil and N_PDS_ParMap are registered to
//                  minimize input for matrix and vector creation
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
// Revision Number: $Revision: 1.69.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:45 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#include <iostream>

#ifdef HAVE_ALGORITHM
#include <algorithm>
#else
#ifdef HAVE_ALGO_H
#include <algo.h>
#else
#error Must have either <algorithm> or <algo.h>!
#endif
#endif
//#if defined LINUX
//#include <algo.h>
//#elif KAILINUX
//#include <algorithm>
//#endif

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>

#include <N_LAS_System.h>

#include <N_TOP_Topology.h>
#include <N_ANP_AnalysisInterface.h>
#include <N_MPDE_Manager.h>

#include <N_PDS_GlobalAccessor.h>
#include <N_PDS_ParMap.h>

#include <N_LAS_QueryUtil.h>
#include <N_LAS_Builder.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>

#include <Epetra_MapColoring.h>

#ifdef Xyce_VERBOSE_LINEAR
 #include <N_ERH_ErrorMgr.h>
#endif

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::~N_LAS_System
// Purpose       : destructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
N_LAS_System::~N_LAS_System()
{
  if ( jacobianMatrixPtr_ ) delete jacobianMatrixPtr_;
  if ( jdxpVectorPtr_ ) delete jdxpVectorPtr_;
  if ( rhsVectorPtr_ ) delete rhsVectorPtr_;
  if ( fVectorPtr_ ) delete fVectorPtr_;
  if ( solnColoringPtr_ ) delete solnColoringPtr_;
  if ( initialConditionColoringPtr_ ) delete initialConditionColoringPtr_;
#ifdef Xyce_DEBUG_VOLTLIM
  if ( dxVoltlimVectorPtr_ ) delete dxVoltlimVectorPtr_;
  // old DAE:
  if ( jdx2VectorPtr_ ) delete jdx2VectorPtr_;
  if ( jacTestMatrixPtr_ ) delete jacTestMatrixPtr_;
  // new DAE:
  if ( Fdx2VectorPtr_ ) delete Fdx2VectorPtr_;
  if ( dFdxTestMatrixPtr_ ) delete dFdxTestMatrixPtr_;
  if ( Qdx2VectorPtr_ ) delete Qdx2VectorPtr_;
  if ( dQdxTestMatrixPtr_ ) delete dQdxTestMatrixPtr_;
#endif
  if (dFdxdVpVectorPtr_) delete dFdxdVpVectorPtr_;
  if (dQdxdVpVectorPtr_) delete dQdxdVpVectorPtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::getGlobalSolutionSize
// Purpose       : returns size of global solution vector
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2/14/06
//-----------------------------------------------------------------------------
int N_LAS_System::getGlobalSolutionSize()
{
  N_LAS_Vector * tmpVec = lasBuilder_->createVector();
  int size = tmpVec->globalLength();
  delete tmpVec;
  return size;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::getGlobalStateSize
// Purpose       : returns size of global solution vector
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/5/10
//-----------------------------------------------------------------------------
int N_LAS_System::getGlobalStateSize()
{
  N_LAS_Vector * tmpVec = lasBuilder_->createStateVector();
  int size = tmpVec->globalLength();
  delete tmpVec;
  return size;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::getSolutionSize
// Purpose       : returns size of local solution vector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
int N_LAS_System::getSolutionSize()
{
  N_LAS_Vector * tmpVec = lasBuilder_->createVector();
  int size = tmpVec->localLength();
  delete tmpVec;
  return size;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::getStateSize
// Purpose       : returns size of local state vector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
int N_LAS_System::getStateSize()
{
  N_LAS_Vector * tmpVec = lasBuilder_->createStateVector();
  int size = tmpVec->localLength();
  delete tmpVec;
  return size;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::getRHSSize
// Purpose       : returns size of local RHS vector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
int N_LAS_System::getRHSSize()
{
  N_LAS_Vector * tmpVec = lasBuilder_->createVector();
  int size = tmpVec->localLength();
  delete tmpVec;
  return size;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::initializeSystem
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
bool N_LAS_System::initializeSystem()
{
  bool bSuccess = true;

//#ifndef Xyce_MPDE
  // ok, this is ugly, but we need to query the MPDE manager to see if we're
  // running in MPDE mode -- in which case coloring doesn't work
  if (anpIntPtr_->getBlockAnalysisFlag() == false) {
//  if(  pdsMgr_->getTopology()->getTIAManager()->getTIAControlAlgorithm()->getMPDEManager()->runFlag() == false ) {
    registerSolnColoring( lasBuilder_->createSolnColoring() );
    registerICColoring( lasBuilder_->createInitialConditionColoring() );
  }
//#endif // Xyce_MPDE

  bSuccess = bSuccess && registerRHSVector( lasBuilder_->createVector() );
  bSuccess = bSuccess && registerJDXPVector( lasBuilder_->createVector() );
  bSuccess = bSuccess && registerFVector( lasBuilder_->createVector() );
#ifdef Xyce_DEBUG_VOLTLIM
  bSuccess = bSuccess && registerDxVoltlimVector( lasBuilder_->createVector() );
 
  // old DAE:
  bSuccess = bSuccess && registerJDX2Vector( lasBuilder_->createVector() );
  bSuccess = bSuccess && registerJacTestMatrix( lasBuilder_->createMatrix() );

  // new DAE:
  bSuccess = bSuccess && registerFDX2Vector( lasBuilder_->createVector() );
  bSuccess = bSuccess && registerdFdxTestMatrix( lasBuilder_->createMatrix() );
  bSuccess = bSuccess && registerQDX2Vector( lasBuilder_->createVector() );
  bSuccess = bSuccess && registerdQdxTestMatrix( lasBuilder_->createMatrix() );
#endif

  // these are needed for new-DAE:
  registerdFdxdVpVector ( lasBuilder_->createVector() );
  registerdQdxdVpVector ( lasBuilder_->createVector() );

  return ( bSuccess && registerJacobianMatrix( lasBuilder_->createMatrix() ) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::updateExternValsSolnVector
// Purpose       : updates off proc values of soln vectors
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/18/00
//-----------------------------------------------------------------------------
bool N_LAS_System::updateExternValsSolnVector(N_LAS_MultiVector * solnVector)
{
  N_PDS_GlobalAccessor * Accessor = pdsMgr_->getGlobalAccessor( "SOLUTION" );

  Accessor->migrateMultiVector(solnVector);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::updateExternValsStateVector
// Purpose       : updates off proc values of state vectors
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/18/00
//-----------------------------------------------------------------------------
bool N_LAS_System::updateExternValsStateVector(N_LAS_MultiVector * stateVector)
{
  N_PDS_GlobalAccessor * Accessor = pdsMgr_->getGlobalAccessor( "STATE" );

  Accessor->migrateMultiVector(stateVector);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::updateExternValsStoreVector
// Purpose       : updates off proc values of store vectors
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
bool N_LAS_System::updateExternValsStoreVector(N_LAS_MultiVector * storeVector)
{
  N_PDS_GlobalAccessor * Accessor = pdsMgr_->getGlobalAccessor( "STORE" );

  Accessor->migrateMultiVector(storeVector);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::localInMatrix
// Purpose       : tests if row is local for jacobian
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/18/00
//-----------------------------------------------------------------------------
bool N_LAS_System::localInMatrix(const int & row) const
{
  return localInSolnVector( row );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::localInSolnVector
// Purpose       : tests if row is local for Soln Vector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/18/00
//-----------------------------------------------------------------------------
bool N_LAS_System::localInSolnVector(const int & row) const
{
#ifdef BAD_STL
	vector<int>::const_iterator iterB = lasQueryUtil_->rowList_GID().begin();
	vector<int>::const_iterator iterE = lasQueryUtil_->rowList_GID().end();
	int total = 0;
	for( ; iterB != iterE ; ++iterB )
  {
		if( *iterB == row )
    {
			++total;
    }
  }
	return total;
#else
  return count(lasQueryUtil_->rowList_GID().begin(),
               lasQueryUtil_->rowList_GID().end(),
               row );
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::localInStateVector
// Purpose       : tests if row is local for State Vector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/18/00
//-----------------------------------------------------------------------------
bool N_LAS_System::localInStateVector(const int & row) const
{
#ifdef BAD_STL
	vector<int>::const_iterator iterB = lasQueryUtil_->rowList_StateGID().begin();
	vector<int>::const_iterator iterE = lasQueryUtil_->rowList_StateGID().end();
	int total = 0;
	for( ; iterB != iterE ; ++iterB )
  {
		if( *iterB == row )
    {
			++total;
    }
  }
	return total;
#else
  return count(lasQueryUtil_->rowList_StateGID().begin(),
               lasQueryUtil_->rowList_StateGID().end(),
               row );
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_System::debug
// Purpose       : debug output
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 11/8/00
//-----------------------------------------------------------------------------
void N_LAS_System::debug() const
{
  cout << "Linear System Debug Output" << endl;
  cout << "--------------------------" << endl;

  cout << "RHS Vector:" << endl;
  rhsVectorPtr_->printPetraObject();

  cout << "Current Solution Vector:" << endl;
  (*currSolVectorPtrPtr_)->printPetraObject();

  cout << "Next Solution Vector:" << endl;
  (*nextSolVectorPtrPtr_)->printPetraObject();

  cout << "Next Solution Derivative Vector:" << endl;
  (*nextSolDerivVectorPtrPtr_)->printPetraObject();

  cout << "Current State Vector:" << endl;
  (*currStaVectorPtrPtr_)->printPetraObject();

  cout << "Next State Vector:" << endl;
  (*nextStaVectorPtrPtr_)->printPetraObject();

  cout << "Temp Solution Vector:" << endl;
  (*tmpSolVectorPtrPtr_)->printPetraObject();

  cout << "Temp Solution Derivative Vector:" << endl;
  (*tmpSolDerivVectorPtrPtr_)->printPetraObject();

  cout << "Temp State Vector:" << endl;
  (*tmpStaVectorPtrPtr_)->printPetraObject();

  cout << "Temp State Derivative Vector:" << endl;
  (*tmpStaDerivVectorPtrPtr_)->printPetraObject();

}

