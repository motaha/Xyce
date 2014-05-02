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
// Filename      : $RCSfile: N_LOA_HBLoader.C,v $
// Purpose       :
// Special Notes :
// Creator       :
// Creation Date :
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.21 $
// Revision Date  : $Date: 2014/02/24 23:49:23 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------

#include <N_LOA_HBLoader.h>
#include <N_MPDE_Discretization.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Builder.h>
#include <N_LAS_HBBuilder.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_Matrix.h>

//#include <N_LAS_Vector.h>
#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_Misc.h>
#include <N_UTL_fwd.h>

#include <N_PDS_ParMap.h>
#include <N_DEV_DeviceInterface.h>

#include <Epetra_MultiVector.h>
#include <Epetra_BlockMap.h>

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerHBBuilder
// Purpose       : Registration method for the HB builder
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, Sandia Labs
// Creation Date : 11/08/13
//-----------------------------------------------------------------------------
void N_LOA_HBLoader::registerHBBuilder( Teuchos::RCP<N_LAS_HBBuilder> hbBuilderPtr )
{
  hbBuilderPtr_ = hbBuilderPtr; 

  // Now initialize all the frequency domain working vectors.
  bXtPtr_ = hbBuilderPtr_->createTimeDomainBlockVector();
  bVtPtr_ = hbBuilderPtr_->createTimeDomainBlockVector();
  bmdQdxPtr_ = rcp_dynamic_cast<N_LAS_BlockMatrix>(rcp(hbBuilderPtr_->createMatrix()));
  bmdFdxPtr_ = rcp_dynamic_cast<N_LAS_BlockMatrix>(rcp(hbBuilderPtr_->createMatrix()));

  bStoreVecFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeStoreBlockVector();
  bStoreLeadCurrQCompVecFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeStoreBlockVector(); 
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::registerAppBuilder
// Purpose       : Registration method for the application builder
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, Sandia Labs
// Creation Date : 11/08/13
//-----------------------------------------------------------------------------
void N_LOA_HBLoader::registerAppBuilder( Teuchos::RCP<N_LAS_Builder> appBuilderPtr )
{
  appBuilderPtr_ = appBuilderPtr;

  // Now initialize all the time domain working vectors.
  appVecPtr_ = rcp(appBuilderPtr_->createVector());
  appNextStaVecPtr_ = rcp(appBuilderPtr_->createStateVector());
  appCurrStaVecPtr_ = rcp(appBuilderPtr_->createStateVector());
  appLastStaVecPtr_ = rcp(appBuilderPtr_->createStateVector());

  appdQdxPtr_ = rcp(appBuilderPtr_->createMatrix());
  appdFdxPtr_ = rcp(appBuilderPtr_->createMatrix());

  appNextStoVecPtr_ = rcp(appBuilderPtr_->createStoreVector());
  appCurrStoVecPtr_ = rcp(appBuilderPtr_->createStoreVector());
  appLastStoVecPtr_ = rcp(appBuilderPtr_->createStoreVector());
  appStoLeadCurrQCompVecPtr_ = rcp(appBuilderPtr_->createStoreVector()); 
}


//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::setFastTimes
// Purpose       : Assign times for fast time scale
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_LOA_HBLoader::setFastTimes( const std::vector<double> & times )
{
  times_ = times;
  constructPeriodicTimes();
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::constructPeriodicTimes
// Purpose       : Make a copy of times_ vector with padding
//                 at the beginning and end to make the calculation
//                 of derivatives easier around the periodicy condition
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_LOA_HBLoader::constructPeriodicTimes()
{
  // we will pad our array of times with 2 times the width, one
  // at the top and one at the bottom of the array
  periodicTimesOffset_ = fastTimeDiscPtr_->Width();
  int timesSize = times_.size();
  period_ = times_[timesSize - 1];
  periodicTimes_.resize( timesSize + 2*periodicTimesOffset_ );
  for( int i=0; i< periodicTimesOffset_; i++ )
  {
    periodicTimes_[i] = times_[i + timesSize - periodicTimesOffset_ - 1] - period_;
  }
  for( int i=periodicTimesOffset_; i< (timesSize + periodicTimesOffset_); i++ )
  {
    periodicTimes_[i] = times_[i - periodicTimesOffset_];
  }
  for( int i=(timesSize+periodicTimesOffset_); i< (timesSize + 2*periodicTimesOffset_); i++ )
  {
    periodicTimes_[i] = times_[i - timesSize - 1] + period_;
  }

}
//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool N_LOA_HBLoader::loadDAEMatrices( N_LAS_Vector * X,
                                     N_LAS_Vector * S,
                                     N_LAS_Vector * dSdt,
                                     N_LAS_Vector * Store,
                                     N_LAS_Matrix * dQdx,
                                     N_LAS_Matrix * dFdx)
{
#ifdef Xyce_DEBUG_HB
#endif // Xyce_DEBUG_HB
  if ( matrixFreeFlag_ )
  {
#ifdef Xyce_DEBUG_HB
    Xyce::dout() << std::endl
           << Xyce::section_divider << std::endl
           << "  N_LOA_HBLoader::loadDAEMatrices:  matrixFree case" << std::endl;
#endif // Xyce_DEBUG_HB
    dQdx->put(0.0);
    dFdx->put(0.0);
    // Do nothing in the Matrix Free Case.
    return(true);
  }
  else
  {
#ifdef Xyce_DEBUG_HB
    Xyce::dout() << std::endl
           << Xyce::section_divider << std::endl
           << "  N_LOA_HBLoader::loadDAEMatrices: Time dependent with matrix case" << std::endl;
#endif // Xyce_DEBUG_HB

    return(loadTimeDepDAEMatrices(X,S,dSdt, Store, dQdx,dFdx));
  }

}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::loadTimeDepDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/03/08
//-----------------------------------------------------------------------------
bool N_LOA_HBLoader::loadTimeDepDAEMatrices( N_LAS_Vector * X,
                                     N_LAS_Vector * S,
                                     N_LAS_Vector * dSdt,
                                     N_LAS_Vector * Store,
                                     N_LAS_Matrix * dQdx,
                                     N_LAS_Matrix * dFdx)
{
#ifdef Xyce_DEBUG_HB
  Xyce::dout() << std::endl
         << Xyce::section_divider << std::endl
         << "  N_LOA_HBLoader::loadTimeDepDAEMatrices" << std::endl;
#endif // Xyce_DEBUG_HB
  //Zero out matrices
  dQdx->put(0.0);
  dFdx->put(0.0);

  N_LAS_Vector appdSdt( *appNextStaVecPtr_ );

  N_LAS_BlockMatrix & bdQdx = *dynamic_cast<N_LAS_BlockMatrix*>(dQdx);
  N_LAS_BlockMatrix & bdFdx = *dynamic_cast<N_LAS_BlockMatrix*>(dFdx);
  N_LAS_BlockVector & bX = *dynamic_cast<N_LAS_BlockVector*>(X);
#ifdef Xyce_FLEXIBLE_DAE_LOADS
  N_LAS_BlockVector & bS = *dynamic_cast<N_LAS_BlockVector*>(S);
  N_LAS_BlockVector & bdSdt = *dynamic_cast<N_LAS_BlockVector*>(dSdt);
  N_LAS_BlockVector & bStore = *dynamic_cast<N_LAS_BlockVector*>(Store);
#endif // Xyce_FLEXIBLE_DAE_LOADS

  int BlockCount = bX.blockCount();
  for( int i = 0; i < BlockCount; ++i )
  {
#ifdef Xyce_DEBUG_HB
    Xyce::dout() << "Processing diagonal matrix block " << i << " of " << BlockCount-1 << std::endl;
#endif // Xyce_DEBUG_HB
#ifdef Xyce_FLEXIBLE_DAE_LOADS
    //Set Time for fast time scale somewhere
    state_.fastTime = times_[i];
    devInterfacePtr_->setFastTime( times_[i] );

    //Update the sources
    appLoaderPtr_->updateSources();

    *appVecPtr_ = bX.block(i);
    *appNextStaVecPtr_ = bS.block(i);
    appdSdt = bdSdt.block(i);
    *appNextStoVecPtr_ = bStore.block(i);

    appLoaderPtr_->loadDAEMatrices( &*appVecPtr_, &*appNextStaVecPtr_, &appdSdt, &*appNextStoVecPtr_, &*appdQdxPtr_, &*appdFdxPtr_);

    bdQdx.block(i,i).add( *appdQdxPtr_ );
    bdFdx.block(i,i).add( *appdFdxPtr_ );
#else
    //For now, the matrices are loaded during the loadDAEVectors method
    //Just copied here
    bdQdx.block(i,i).add( bmdQdxPtr_->block(i,i) );
    bdFdx.block(i,i).add( bmdFdxPtr_->block(i,i) );

#endif // Xyce_FLEXIBLE_DAE_LOADS
  }

  // Fast Time scale discretization terms:
  // These are d(dQ/dt1)/dx terms, but go into bdFdx.
  // For this procedure, need to re-use the app matrix, appdQdx.
  N_LAS_Matrix & dQdxdt = *appdQdxPtr_;

  int Start = fastTimeDiscPtr_->Start();
  int Width = fastTimeDiscPtr_->Width();

  const std::vector<double> & Coeffs = fastTimeDiscPtr_->Coeffs();

  for( int i = 0; i < BlockCount; ++i )
  {
#ifdef Xyce_DEBUG_HB
    Xyce::dout() << "Processing off diagonal matrix blocks on row " << i << " of " << BlockCount-1 << std::endl;
#endif // Xyce_DEBUG_HB
    int Loc;
    int indexT1 = i + Start + periodicTimesOffset_;
    int indexT2 = indexT1 + Width - 1;
    double invh2 = 1.0 / (periodicTimes_[indexT2] - periodicTimes_[indexT1]);

    for( int j = 0; j < Width; ++j )
    {
      Loc = i + (j+Start);

      if( Loc < 0 )
      {
        Loc += BlockCount;
      }
      else if( Loc > (BlockCount-1) )
      {
        Loc -= BlockCount;
      }

      dQdxdt.put(0.0);
      dQdxdt.add( bdQdx.block(Loc,Loc) );
      dQdxdt.scale( Coeffs[j]*invh2 );
      bdFdx.block(i,Loc).add( dQdxdt );
    }
  }

  dQdx->fillComplete();
  dFdx->fillComplete();

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "HB bX:" << std::endl;
  bX.printPetraObject(std::cout);
  Xyce::dout() << "HB bdQdx:" << std::endl;
  bdQdx.printPetraObject(std::cout);
  Xyce::dout() << "HB bdFdx:" << std::endl;
  bdFdx.printPetraObject(std::cout);
#ifdef Xyce_FLEXIBLE_DAE_LOADS
  Xyce::dout() << "HB bS:" << std::endl;
  bS.printPetraObject(std::cout);
  Xyce::dout() << "HB dSdt:" << std::endl;
  bdSdt.printPetraObject(std::cout);
  Xyce::dout() << "HB bStore:" << std::endl;
  bStore.printPetraObject(std::cout);
#endif // Xyce_FLEXIBLE_DAE_LOADS

  Xyce::dout() << Xyce::section_divider << std::endl;
#endif // Xyce_DEBUG_HB
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::applyDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool N_LOA_HBLoader::applyDAEMatrices( N_LAS_Vector * Xf,
                                     N_LAS_Vector * S,
                                     N_LAS_Vector * dSdt,
                                     N_LAS_Vector * Store,
                                     const N_LAS_Vector & Vf,
                                     N_LAS_Vector * dQdxV,
                                     N_LAS_Vector * dFdxV )
{
  if ( !matrixFreeFlag_ )
  {
    std::string msg="N_LOA_HBLoader::applyDAEMatrices.  This function should only be called in the matrix free case.";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }
#ifdef Xyce_DEBUG_HB
  Xyce::dout() << std::endl
         << Xyce::section_divider << std::endl
         << "  N_LOA_HBLoader::applyDAEMatrices" << std::endl;
#endif // Xyce_DEBUG_HB
  //Zero out matvec vectors
  dQdxV->putScalar(0.0); // This is dQdx * V on output (only used in HB-env simulation)
  dFdxV->putScalar(0.0); // This is dFdx * V on output (This is main output for HB)

  bXtPtr_->putScalar(0.0);
  bVtPtr_->putScalar(0.0);

  N_LAS_BlockVector & bXf = *dynamic_cast<N_LAS_BlockVector*>(Xf);
  // We have to do something special with Vf because AztecOO (or Belos)
  // probably used the Epetra_LinearProblem's Epetra_Maps to create the input
  // vector here.  In this case, Vf is just an N_LAS_Vector and not an
  // N_LAS_BlockVector.
  const N_LAS_BlockVector bVf(Vf, bXf.blockSize());

  permutedIFT(bXf, &*bXtPtr_);
  permutedIFT(bVf, &*bVtPtr_);

  N_LAS_BlockVector & bX = *bXtPtr_;

  N_LAS_BlockVector * bdQdxV = dynamic_cast<N_LAS_BlockVector*>(dQdxV);
  N_LAS_BlockVector * bdFdxV = dynamic_cast<N_LAS_BlockVector*>(dFdxV);

  Teuchos::RCP<N_LAS_BlockVector>  bdQdxVt = hbBuilderPtr_->createTimeDomainBlockVector();
  Teuchos::RCP<N_LAS_BlockVector>  bdFdxVt = hbBuilderPtr_->createTimeDomainBlockVector();

  int BlockCount = bX.blockCount();
  for( int i = 0; i < BlockCount; ++i )
  {
#ifdef Xyce_DEBUG_HB
    Xyce::dout() << "Processing diagonal matrix block " << i << " of " << BlockCount-1 << std::endl;
#endif // Xyce_DEBUG_HB

    // Get the already stored time-domain Jacobian matrices
    appdQdxPtr_ = vecAppdQdxPtr_[i];
    appdFdxPtr_ = vecAppdFdxPtr_[i];

#ifdef Xyce_DEBUG_HB
    Xyce::dout() << "bVtPtr block i = " << i << " : " << std::endl;
    bVtPtr_->block(i).printPetraObject(std::cout);

    Xyce::dout() << "appdQdxPtr_ = " << i << " : " <<  std::endl;
    appdQdxPtr_->printPetraObject(std::cout);

    Xyce::dout() << "appdFdxPtr_ = " << i << " : " <<  std::endl;
    appdFdxPtr_->printPetraObject(std::cout);
#endif // Xyce_DEBUG_HB

    appdQdxPtr_->matvec(false, bVtPtr_->block(i), bdQdxVt->block(i));
    appdFdxPtr_->matvec(false, bVtPtr_->block(i), bdFdxVt->block(i));

#ifdef Xyce_DEBUG_HB
    Xyce::dout() << "bdQdxVt block i = " << i << " : " << std::endl;
    bdQdxVt->block(i).printPetraObject(std::cout);

    Xyce::dout() << "bdFdxVt block i = " << i << " : " << std::endl;
    bdFdxVt->block(i).printPetraObject(std::cout);
#endif // Xyce_DEBUG_HB
  }

  permutedFFT(*bdQdxVt, bdQdxV);
  permutedFFT(*bdFdxVt, bdFdxV);

  int blockCount = bXf.blockCount();
  int blockSize = bXf.block(0).globalLength();
  double omega = 2.0 * M_PI/ period_;

  for( int i = 0; i < blockCount; ++i )
  {
    N_LAS_Vector QVec(bXf.block(i));
    N_LAS_Vector freqVec = bdQdxV->block(i);

    // Only one processor owns each block of the frequency-domain vector
    if (freqVec.localLength() > 0)
    {
      QVec[0] = -freqVec[1]*0.0*omega;
      QVec[1] = freqVec[0]*0.0*omega;

      for (int j=1; j < (blockSize/2+1)/2; ++j)
      {
        QVec[2*j] = -freqVec[2*j+1]*j*omega;
        QVec[2*(blockSize/2-j)] = -freqVec[2*j+1]*j*omega;

        QVec[2*j+1] = freqVec[2*j]*j*omega;
        QVec[2*(blockSize/2-j)+1] = -freqVec[2*j]*j*omega;
      }
    }

    bdFdxV->block(i).update(1.0, QVec , 1.0);
  }

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "HB bX:" << std::endl;
  bX.printPetraObject(std::cout);
  Xyce::dout() << "HB bdQdxV:" << std::endl;
  bdQdxV->printPetraObject(std::cout);
  Xyce::dout() << "HB bdFdxV:" << std::endl;
  bdFdxV->printPetraObject(std::cout);

  Xyce::dout() << Xyce::section_divider << std::endl;
#endif // Xyce_DEBUG_HB
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::updateState
// Purpose       :
// Special Notes : ERK.  This function needs to be a no-op.  The reason
//                 is that the state information needs to be the same
//                 at the time of updateState, loadDAEVectors and
//                 loadDAEMatrices.  Thus, they have to all happen inside
//                 of the same "fast time" loop.  So, this functionality
//                 has been moved up into loadDAEVectors.
//
//                 Note: for similar reasons, loadDAEMatrices is called from
//                 within that function as well.
//
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
bool N_LOA_HBLoader::updateState
 (N_LAS_Vector * nextSolVectorPtr,
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
  bool bsuccess = true;

  // For HB case, this needs to be a no-op.

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool N_LOA_HBLoader::loadDAEVectors( N_LAS_Vector * Xf,
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
                                  N_LAS_Vector * dQdxdVp )
{

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << std::endl
         << Xyce::section_divider << std::endl
         << "  N_LOA_HBLoader::loadDAEVectors" << std::endl;
#endif // Xyce_DEBUG_HB

  //Zero out vectors
  appVecPtr_->putScalar(0.0);
  appNextStaVecPtr_->putScalar(0.0);
  appCurrStaVecPtr_->putScalar(0.0);
  appLastStaVecPtr_->putScalar(0.0);
  N_LAS_Vector appdSdt( *appNextStaVecPtr_ );

  appNextStoVecPtr_->putScalar(0.0);
  appCurrStoVecPtr_->putScalar(0.0);
  appLastStoVecPtr_->putScalar(0.0);
  appStoLeadCurrQCompVecPtr_->putScalar(0.0);

  N_LAS_Vector appQ( *appVecPtr_ );
  N_LAS_Vector appF( *appVecPtr_ );

  N_LAS_Vector appdFdxdVp( *appVecPtr_ );
  N_LAS_Vector appdQdxdVp( *appVecPtr_ );

  // This is a temporary load storage vector.
  N_LAS_Vector dQdt2( *appVecPtr_ );

  bXtPtr_->putScalar(0.0);

  N_LAS_BlockVector & bXf = *dynamic_cast<N_LAS_BlockVector*>(Xf);

  permutedIFT(bXf, &*bXtPtr_);

  // 12/8/06 tscoffe:   Note:  "b" at beginning of variable name means N_LAS_BlockVector
  N_LAS_BlockVector & bX = *bXtPtr_;
  N_LAS_BlockVector & bS = *dynamic_cast<N_LAS_BlockVector*>(S);
  N_LAS_BlockVector & bcurrS = *dynamic_cast<N_LAS_BlockVector*>(currS);
  N_LAS_BlockVector & blastS = *dynamic_cast<N_LAS_BlockVector*>(lastS);
  N_LAS_BlockVector & bdSdt = *dynamic_cast<N_LAS_BlockVector*>(dSdt);
  N_LAS_BlockVector & bStore = *dynamic_cast<N_LAS_BlockVector*>(Store);
  N_LAS_BlockVector & bcurrStore = *dynamic_cast<N_LAS_BlockVector*>(currStore);
  N_LAS_BlockVector & blastStore = *dynamic_cast<N_LAS_BlockVector*>(lastStore);
  N_LAS_BlockVector * bQ = dynamic_cast<N_LAS_BlockVector*>(Q);
  N_LAS_BlockVector * bF = dynamic_cast<N_LAS_BlockVector*>(F);

  N_LAS_BlockVector * bdFdxdVp = dynamic_cast<N_LAS_BlockVector*>(dFdxdVp);
  N_LAS_BlockVector * bdQdxdVp = dynamic_cast<N_LAS_BlockVector*>(dQdxdVp);


//  N_LAS_Vector * storeLeadCurrQComp,

  N_LAS_BlockVector & bstoreLeadCurrQComp = *dynamic_cast<N_LAS_BlockVector*>(storeLeadCurrQComp);

  N_LAS_BlockVector  bQt(*bXtPtr_);
  N_LAS_BlockVector  bFt(*bXtPtr_);

  N_LAS_BlockVector  bdQdxdVpt(*bXtPtr_);
  N_LAS_BlockVector  bdFdxdVpt(*bXtPtr_);


//  N_LAS_BlockVector  bstoreLeadCurrQCompt(*bXtPt  );

#ifndef Xyce_FLEXIBLE_DAE_LOADS
  bmdQdxPtr_->put(0.0);
  bmdFdxPtr_->put(0.0);
#endif

  int BlockCount = bX.blockCount();

  // We are storing the time domain Jacobians, initialize the memory here 
  if ((int)vecAppdQdxPtr_.size()!=BlockCount)
  {
    vecAppdQdxPtr_.resize( BlockCount );
    vecAppdFdxPtr_.resize( BlockCount );
    for ( int i=0; i < BlockCount; ++i )
    {
      vecAppdQdxPtr_[i] = rcp(appBuilderPtr_->createMatrix());
      vecAppdFdxPtr_[i] = rcp(appBuilderPtr_->createMatrix());
    }
  }

  // Now perform implicit application of frequency domain Jacobian. 
  for( int i = 0; i < BlockCount; ++i )
  {
#ifdef Xyce_DEBUG_HB
    Xyce::dout() << "Processing vectors for block " << i << " of " << BlockCount-1 << std::endl;
#endif // Xyce_DEBUG_HB
    //Set Time for fast time scale somewhere
    state_.fastTime = times_[i];
    devInterfacePtr_->setFastTime( times_[i] );

#ifdef Xyce_DEBUG_HB
    Xyce::dout() << "Calling updateSources on the appLoader" << std::endl;
#endif // Xyce_DEBUG_HB
    //Update the sources
    appLoaderPtr_->updateSources();  // this is here to handle "fast" sources.

    *appVecPtr_ = bX.block(i);
    *appNextStaVecPtr_ = bS.block(i);
    *appCurrStaVecPtr_ = bcurrS.block(i);
    *appLastStaVecPtr_ = blastS.block(i);
    appdSdt = bdSdt.block(i);
    *appNextStoVecPtr_ = bStore.block(i);
    *appCurrStoVecPtr_ = bcurrStore.block(i);
    *appLastStoVecPtr_ = blastStore.block(i);
    *appStoLeadCurrQCompVecPtr_ = bstoreLeadCurrQComp.block(i);

//blastStore.block(i); // Need to update to correct block!

#ifdef Xyce_DEBUG_HB
    Xyce::dout() << "Updating State for block " << i << " of " << BlockCount-1 << std::endl;
#endif // Xyce_DEBUG_HB

    // Note: This updateState call is here (instead of in the
    // N_LOA_HBLoader::updateState function) because it has to be called
    // for the same fast time point.
    appLoaderPtr_->updateState
      ( &*appVecPtr_,
        &*appVecPtr_,  // note, this is a placeholder! ERK
        &*appVecPtr_,  // note, this is a placeholder! ERK
        &*appNextStaVecPtr_, &*appCurrStaVecPtr_ , &*appLastStaVecPtr_,
        &*appNextStoVecPtr_, &*appCurrStoVecPtr_ , &*appLastStoVecPtr_
        );

    bS.block(i) = *appNextStaVecPtr_;
    bcurrS.block(i) = *appCurrStaVecPtr_;
    blastS.block(i) = *appLastStaVecPtr_;
    bStore.block(i) = *appNextStoVecPtr_;
    bcurrStore.block(i) = *appCurrStoVecPtr_;
    blastStore.block(i) = *appLastStoVecPtr_;

#ifdef Xyce_DEBUG_HB
    Xyce::dout() << "Calling loadDAEVectors on the appLoader" << std::endl;
#endif // Xyce_DEBUG_HB

    // This has to be done because the app loader does NOT zero these vectors out.
    appQ.putScalar(0.0);
    appF.putScalar(0.0);
    appdFdxdVp.putScalar(0.0);
    appdQdxdVp.putScalar(0.0);

    appLoaderPtr_->loadDAEVectors
      ( &*appVecPtr_,
        &*appVecPtr_,  // note, this is a placeholder! ERK
        &*appVecPtr_,  // note, this is a placeholder! ERK
        &*appNextStaVecPtr_, &*appCurrStaVecPtr_, &*appLastStaVecPtr_, &appdSdt,
        &*appNextStoVecPtr_, &*appCurrStoVecPtr_, &*appLastStoVecPtr_, &*appStoLeadCurrQCompVecPtr_,
        &appQ, &appF, &appdFdxdVp, &appdQdxdVp );

    bQt.block(i) = appQ;
    bFt.block(i) = appF;

    bdQdxdVpt.block(i) = appdQdxdVp;
    bdFdxdVpt.block(i) = appdFdxdVp;

    bstoreLeadCurrQComp.block(i) = *appStoLeadCurrQCompVecPtr_;

    bStore.block(i) = *appNextStoVecPtr_;  //lead current get loaded during loadDAEVectors

//    Xyce::dout() << "HB Loader Store Vector TD" << std::endl;
//    bStore.printPetraObject(std::cout);

//    Xyce::dout() << "HB Loader Store  Q Vector TD" << std::endl;
//    bstoreLeadCurrQComp.printPetraObject(std::cout);

    // Store the time domain Jacobian for future use.
    vecAppdQdxPtr_[i]->put(0.0);
    vecAppdFdxPtr_[i]->put(0.0);

    // Load dQdx and dFdx into the storage location for this time point.
    appLoaderPtr_->loadDAEMatrices( &*appVecPtr_, &*appNextStaVecPtr_, &appdSdt, &*appNextStoVecPtr_, &*vecAppdQdxPtr_[i],  &*vecAppdFdxPtr_[i]);
  }

  permutedFFT(bQt, bQ);

  permutedFFT(bFt, bF);

  permutedFFT(bdQdxdVpt, bdQdxdVp);

  permutedFFT(bdFdxdVpt, bdFdxdVp);

  bStoreVecFreqPtr_->putScalar(0.0);
  bStoreLeadCurrQCompVecFreqPtr_->putScalar(0.0);

  permutedFFT(bStore, &*bStoreVecFreqPtr_);
  permutedFFT(bstoreLeadCurrQComp, &*bStoreLeadCurrQCompVecFreqPtr_);  

//  Xyce::dout() << "HB Store Vector FD" << std::endl;
//  bStoreVecFreqPtr_->printPetraObject(std::cout);
//  bStoreLeadCurrQCompVecFreqPtr_->printPetraObject(std::cout);

  int blockCount = bXf.blockCount();
  int blockSize = bXf.block(0).globalLength();
  double omega = 2.0 * M_PI/ period_;

  for( int i = 0; i < blockCount; ++i )
  {
    // Create work vectors from the current frequency block vector
    // NOTE:  This needs to be done for each block to make sure that the
    //        map is the same as the bF block.
    N_LAS_Vector QVec(bQ->block(i));
    N_LAS_Vector freqVec = bQ->block(i);

    N_LAS_Vector dQdxdVpVec(bdQdxdVp->block(i));
    N_LAS_Vector freqVec1 = bdQdxdVp->block(i);

    // Only one processor owns each block of the frequency-domain vector
    if (freqVec.localLength() > 0)
    {
      QVec[0] = -freqVec[1]*0.0*omega;
      QVec[1] = freqVec[0]*0.0*omega;

      dQdxdVpVec[0] = -freqVec1[1]*0.0*omega;
      dQdxdVpVec[1] = freqVec1[0]*0.0*omega;

      for (int j=1; j < (blockSize/2+1)/2; ++j)
      {
        QVec[2*j] = -freqVec[2*j+1]*j*omega;
        QVec[2*(blockSize/2-j)] = -freqVec[2*j+1]*j*omega;

        QVec[2*j+1] = freqVec[2*j]*j*omega;
        QVec[2*(blockSize/2-j)+1] = -freqVec[2*j]*j*omega;

        dQdxdVpVec[2*j] = -freqVec1[2*j+1]*j*omega;
        dQdxdVpVec[2*(blockSize/2-j)] = -freqVec1[2*j+1]*j*omega;

        dQdxdVpVec[2*j+1] = freqVec1[2*j]*j*omega;
        dQdxdVpVec[2*(blockSize/2-j)+1] = -freqVec1[2*j]*j*omega;

      }
    }

    bF->block(i).update(1.0, QVec , 1.0);

    bdFdxdVp->block(i).update(1.0, dQdxdVpVec, 1.0);

  }


  blockCount = bStoreVecFreqPtr_->blockCount();
  blockSize =  bStoreVecFreqPtr_->blockSize();

  for( int i = 0; i < blockCount; ++i )
  {
    // Create work vectors from the current frequency block vector
    // NOTE:  This needs to be done for each block to make sure that the
    //        map is the same as the bF block.

    N_LAS_Vector stoLeadCurrdQdtVec(bStoreLeadCurrQCompVecFreqPtr_->block(i));
    N_LAS_Vector stoLeadCurrQVec = bStoreLeadCurrQCompVecFreqPtr_->block(i);


    // Only one processor owns each block of the frequency-domain vector
    if (stoLeadCurrQVec.localLength() > 0)
    {

      stoLeadCurrdQdtVec[0] = stoLeadCurrQVec[1]*0.0*omega;
      stoLeadCurrdQdtVec[1] = stoLeadCurrQVec[0]*0.0*omega;

      for (int j=1; j < (blockSize/2+1)/2; ++j)
      {

        stoLeadCurrdQdtVec[2*j] = -stoLeadCurrQVec[2*j+1]*j*omega;
        stoLeadCurrdQdtVec[2*(blockSize/2-j)] = -stoLeadCurrQVec[2*j+1]*j*omega;

        stoLeadCurrdQdtVec[2*j+1] = stoLeadCurrQVec[2*j]*j*omega;
        stoLeadCurrdQdtVec[2*(blockSize/2-j)+1] = -stoLeadCurrQVec[2*j]*j*omega;
      }
    }

    bStoreVecFreqPtr_->block(i).update(1.0, stoLeadCurrdQdtVec, 1.0);

//    Xyce::dout() << "HB Store Vector  dqdt + f(x) FD" << std::endl;
//    bStoreVecFreqPtr_->printPetraObject(std:: cout);
  }
//  permutedIFT(*bStoreVecFreqPtr_, &bStore);
#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "HB X Vector" << std::endl;
  bX.printPetraObject(std::cout);
 //  Xyce::dout() << "HB S Vector" << std::endl;
 //  bS.printPetraObject(Xyce::dout());
 //  Xyce::dout() << "HB dSdt Vector" << std::endl;
 //  bdSdt.printPetraObject(Xyce::dout());
  Xyce::dout() << "HB Store Vector" << std::endl;
  bStore.printPetraObject(std::cout);
  Xyce::dout() << "HB Q Vector" << std::endl;
  bQ->printPetraObject(std::cout);
  Xyce::dout() << "HB F Vector" << std::endl;
  bF->printPetraObject(std::cout);
  Xyce::dout() << "HB bdFdxdVp Vector" << std::endl;
  bdFdxdVp->printPetraObject(std::cout);
  Xyce::dout() << "HB bdQdxdVp Vector" << std::endl;
  bdQdxdVp->printPetraObject(std::cout);

#ifndef Xyce_FLEXIBLE_DAE_LOADS
  Xyce::dout() << "HB bmdQdx_" << std::endl;
  bmdQdxPtr_->printPetraObject(std::cout);
  Xyce::dout() << "HB bmdFdx_" << std::endl;
  bmdFdxPtr_->printPetraObject(std::cout);
#endif // Xyce_FLEXIBLE_DAE_LOADS
  Xyce::dout() << Xyce::section_divider << std::endl;
#endif // Xyce_DEBUG_HB

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::permutedFFT
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 09/05/08
//---------------------------------------------------------------------------
void N_LOA_HBLoader::permutedFFT(const N_LAS_BlockVector & xt, N_LAS_BlockVector * xf)
{
  int blockCount = xt.blockCount();
  int N = xt.block(0).globalLength();

  // It's necessary to get the blockmap from Epetra because the N_PDS_ParMap is not always guaranteed to be valid.
  Epetra_BlockMap blockMap = xt.block(0).epetraObj().Map();

  N_UTL_FFTInterface<std::vector<double> > HBTransform(blockCount);

  // Register vectors with the FFT interface.
  std::vector<double> inputSignal(blockCount, 0.0);
  std::vector<double> outputSignal((blockCount+1), 0.0);
  std::vector<double> temp(0);  // Dummy vector, not needed since we only compute FFT.
  HBTransform.registerVectors( inputSignal, &outputSignal, temp, &temp );

  // Loop through all the solution variables and compute the FFT of each variable.
  for (int j=0; j<N; j++)
  {
    // Get the local id for this variable.
    int lid = blockMap.LID(j);

    N_LAS_Vector& freqVecRef = xf->block(j);

#ifdef Xyce_DEBUG_HB
    std::cout << "Xf block(" << j << ") before = " << std::endl;
    freqVecRef.printPetraObject(std::cout);
#endif // Xyce_DEBUG_HB

    for (int i=0; i<blockCount; ++i)
    {
      N_LAS_Vector& timeVecRef = xt.block(i);

      if (lid >= 0)
      {
        inputSignal[i] = timeVecRef[lid];
#ifdef Xyce_DEBUG_HB
        std::cout << "Block " << j << ": inputSignal ("<<i<<") = " << inputSignal[i] << std::endl;
#endif // Xyce_DEBUG_HB
      }
    }

    // Only perform the FFT on the processor that owns the solution variable.
    if (lid >= 0)
    {
      // Calculate the FFT for the inputSignal.
      HBTransform.calculateFFT();

#ifdef Xyce_DEBUG_HB
      for (int i=0; i<blockCount+1; ++i)
      {
        std::cout << "Block " << j << ": outputSignal ("<<i<<") = " << outputSignal[i] << std::endl;
      }
#endif // Xyce_DEBUG_HB

      freqVecRef[0] =  outputSignal[0]/blockCount;
      freqVecRef[1] =  outputSignal[1]/blockCount;

      for (int i=1; i<(blockCount+1)/2; ++i)
      {
        freqVecRef[2*i] =  outputSignal[2*i]/blockCount;
        freqVecRef[2*(blockCount-i)] =  outputSignal[2*i]/blockCount;

        freqVecRef[2*i+1] =  outputSignal[2*i+1]/blockCount;
        freqVecRef[2*(blockCount-i)+1] = -outputSignal[2*i+1]/blockCount;
      }
    }

#ifdef Xyce_DEBUG_HB
    std::cout << "Xf block(" << j << ") after = " << std::endl;
    freqVecRef.printPetraObject(std::cout);
#endif // Xyce_DEBUG_HB

  }

}

//-----------------------------------------------------------------------------
// Function      : N_LOA_HBLoader::permutedIFT
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 09/05/08
//---------------------------------------------------------------------------
void  N_LOA_HBLoader::permutedIFT(const N_LAS_BlockVector & xf, N_LAS_BlockVector * xt)
{
  int blockCount = xf.blockCount();
  int N = xf.block(0).globalLength();

  // It's necessary to get the blockmap from Epetra because the N_PDS_ParMap is not always guaranteed to be valid.
  Epetra_BlockMap blockMap = (xt->block(0)).epetraObj().Map();

#ifdef Xyce_DEBUG_HB
  std::cout << "FD blockCount= " << blockCount   << ", fN =" << N <<  std::endl;

  int tmpblockCount = xt->blockCount();
  int tmpN = xt->block(0).globalLength();
  std::cout << "TD blockCount= " << tmpblockCount   << ", tN =" << tmpN <<  std::endl;
#endif // Xyce_DEBUG_HB

  N_UTL_FFTInterface<std::vector<double> > HBTransform((N/2));

  // Register input and output signals for the IFT.
  std::vector<double> inputSignal((N/2+1), 0.0);
  std::vector<double> outputSignal((N/2), 0.0);
  std::vector<double> temp(0);  // Dummy vector, not needed since we only compute IFT.
  HBTransform.registerVectors( temp, &temp, inputSignal, &outputSignal );

  for (int j=0; j<blockCount; j++)
  {
    N_LAS_Vector& freqVecRef = xf.block(j);

    if (freqVecRef.localLength() > 0)
    {
      for (int i=0; i<(N/2+1); ++i)
      {
        inputSignal[i] = freqVecRef[i];

#ifdef Xyce_DEBUG_HB
        std::cout << "Block " << j << ": inputSignal i= " << i << ", input= " << inputSignal[i] << std::endl;
#endif // Xyce_DEBUG_HB
      }

      // Calculate the inverse FFT.
      HBTransform.calculateIFT();
    }

    for (int i=0; i<(N/2); ++i)
    {
      N_LAS_Vector& timeVecRef = xt->block(i);
      if (freqVecRef.localLength() > 0)
      {
        timeVecRef[blockMap.LID(j)] =  outputSignal[i]*(N/2);
#ifdef Xyce_DEBUG_HB
        std::cout << "Block " << j << ": outputSignal i= " << i << ", output= " << timeVecRef[blockMap.LID(j)] << std::endl;
#endif // Xyce_DEBUG_HB
      }
    }
  }

}

